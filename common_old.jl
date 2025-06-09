import Trixi: flux_hllc
# use the 2D implementation since the 1D version doesn't account for n::Integer < 0
function flux_hllc(u_ll, u_rr, n::SVector{1}, equations::CompressibleEulerEquations1D)
    f = flux_hllc(SVector(u_ll[1], u_ll[2], 0, u_ll[3]), 
                  SVector(u_rr[1], u_rr[2], 0, u_rr[3]), 
                  SVector(n[1], 0.0), 
                  CompressibleEulerEquations2D(equations.gamma))    
    return SVector(f[1], f[2], f[4])
end

function store_fourier_coeffs!(target, param)
    if eltype(param) == Float64
        target .= real.(fft(param))
    else
        target .= SVector{length(eltype(param)), Float64}.([real.(fft(getindex.(param, i))) for i in 1:length(eltype(param))]...)
    end
end

function add_inverse_fourier!(target, param)
    if eltype(param) == Float64
        target .+= real.(ifft(param))
    else
        target .+= SVector{length(eltype(param)), Float64}.([real.(ifft(getindex.(param, i))) for i in 1:length(eltype(param))]...)
    end
end

function rhs!(du, u, cache, t)
    (; weird_Q_skew_nz, M, psi, alpha, dt, blend, entropy_inequality, entropy_blend, blending_strat, filter_strength, volume_flux, low_order_volume_flux, equations, r_H, r_L, r_entropy_rhs, a, θ, v, knapsack_solver!, bc, cube_space, FH_ij_storage, FL_ij_storage, knapsack_shock_capturing) = cache
    fill!(du, zero(eltype(du)))
    fill!(r_H, zero(eltype(du)))
    fill!(r_L, zero(eltype(du)))
    fill!(r_entropy_rhs, zero(eltype(x)))

    if !isnothing(bc)
        u[1] = bc[:, 1]
        u[end] = bc[:, end]
    end

    @. v = cons2entropy.(u, equations)

    for (j, row) in enumerate(weird_Q_skew_nz)
        if entropy_inequality == :semi_local
            if entropy_blend == :grouped || entropy_blend == :viscosity
                a = @view r_entropy_rhs[1:length(row)] # we dont need r_entropy_rhs for this strategy, so store a here
                θ_local = @view θ[1:length(row)] # same for θ
                b = 0.
            elseif entropy_blend == :free
                a = @view r_entropy_rhs[1:length(row) + 1] # we dont need r_entropy_rhs for this strategy, so store a here
                a[end] = zero(eltype(a))
                θ_local = @view θ[1:length(row) + 1] # same for θ
                b = 0.
            end
        end
        for (index, (i, q)) in enumerate(row)
            if abs(q) > 1e-12
                nij = q
                FH_ij_storage[i, j] = volume_flux(u[i], u[j], SVector{1}(nij / norm(nij)), equations)
                
                if entropy_blend == :viscosity
                    FL_ij_storage[i, j] = u[i] - u[j]
                else
                    FL_ij_storage[i, j] = low_order_volume_flux(u[i], u[j], SVector{1}(nij / norm(nij)), equations)
                end

                FH_ij = norm(nij) * FH_ij_storage[i, j]
                FL_ij = norm(nij) * FL_ij_storage[i, j]
                entropy_ij = norm(nij) * (v[j]'FH_ij_storage[i, j] - (psi(u[j]) - psi(u[i])) * nij / norm(nij))
                entropy_lom_ij = norm(nij) * (v[j]'FL_ij_storage[i, j] - (psi(u[j]) - psi(u[i])) * nij / norm(nij))

                if entropy_inequality == :local
                    λ = 0.
                    if abs(v[i]' * (FH_ij - FL_ij)) > 100 * eps() && -v[i]'FH_ij + entropy_ij > 100 * eps()
                        λ = clamp((v[i]'FH_ij - entropy_ij) / (v[i]' * (FH_ij - FL_ij)), 0., 1.)
                    end

                    r_H[i] += λ * FL_ij + (1 - λ) * FH_ij
                    r_H[j] -= λ * FL_ij + (1 - λ) * FH_ij
                elseif entropy_inequality == :semi_local
                    if entropy_blend == :grouped
                        a[index] = v[i]' * (FL_ij - FH_ij) - (entropy_lom_ij - entropy_ij)
                    elseif entropy_blend == :free
                        a[index] = v[i]' * (FL_ij - FH_ij)
                        a[end] += entropy_ij - entropy_lom_ij
                    elseif entropy_blend == :viscosity
                        a[index] = (v[i] - v[j])'FL_ij
                    end

                    if knapsack_solver!.direction == 1
                        # Minimization problem
                        b += -v[i]'FH_ij + entropy_ij
                    else
                        # Maximization problem
                        b += -entropy_lom_ij + v[i]'FL_ij
                    end

                    r_H[i] += FH_ij
                    r_H[j] -= FH_ij
                else
                    r_H[i] += FH_ij
                    r_H[j] -= FH_ij

                    r_entropy_rhs[i] += entropy_ij
                    r_entropy_rhs[j] -= entropy_ij
                end

                r_L[i] += FL_ij
                r_L[j] -= FL_ij
            end
        end

        if entropy_inequality == :semi_local
            if blend == :knapsack
                # Solve a knapsack problem
                knapsack_solver!(θ_local, a, b)

                # Shock capturing only works for maximizers
                if knapsack_solver!.direction == -1 && knapsack_shock_capturing >= 1
                    @. a *= knapsack_shock_capturing * (1 - θ_local)^2 + 1

                    knapsack_solver!(θ_local, a, b)
                end

            elseif blend == :scalar
                # Sort of like elementwise limiting!
                if entropy_blend == :free
                    # Keep the entropy component seperate, do a knapsack problem of size 2
                    tiny_a = @view a[end - 1:end]
                    tiny_θ = @view θ_local[end - 1:end]
                    tiny_a[1] = sum(@view a[1:end - 1])
                    knapsack_solver!(tiny_θ, tiny_a, b)

                    # Now we only need the blending part of tiny_θ
                    θ_local .= tiny_θ[1]
                else
                    sum_a = sum(a)
                    if abs(sum_a) < 100 * eps() || b < 100 * eps() # already satisfies entropy inequality
                        θ_local .= 0.
                    elseif entropy_blend == :viscosity
                        # Unbounded
                        θ_local .= b * a / (a'a)
                    else
                        θ_local .= clamp(b / sum_a, 0., 1.)
                    end
                end
            end

            # Save in corresponding index in cube_space
            for (index, (i, _)) in enumerate(row)
                cube_space[i, j] = θ_local[index]
            end
        end
    end

    if entropy_inequality == :semi_local
        for (j, row) in enumerate(weird_Q_skew_nz)
            for (index, (i, q)) in enumerate(row)
                if i > j
                    # Relax blending coefficients to satisfy semi local entropy inequality for each node
                    θ = max(cube_space[i, j])

                    if alpha >= 0
                        # Let's also satisfy a positivity constraint
                        bound_1 = 0.
                        if r_H[i][1] - r_L[i][1] > 100 * eps()
                            bound_1 = 1 - (1 - alpha) / dt * (M[i, i] * u[i][1] - dt * r_L[i][1]) / (r_H[i][1] - r_L[i][1])
                        end

                        bound_3 = 0.
                        if r_H[i][3] - r_L[i][3] > 100 * eps()
                            bound_3 = 1 - (1 - alpha) / dt * (M[i, i] * u[i][3] - dt * r_L[i][3]) / (r_H[i][3] - r_L[i][3])
                        end
                        
                        θ = max(θ, min(1, bound_1), min(1, bound_3))
                    end

                    cube_space[i, j] = θ
                end
            end
        end

        # Swap pointers for FH_ij_storage and FL_ij_storage pointers if necessary
        if knapsack_solver!.direction == -1
            temp_pointer = FL_ij_storage
            FL_ij_storage = FH_ij_storage
            FH_ij_storage = temp_pointer
        end
        
        if entropy_blend == :viscosity
            for (j, row) in enumerate(weird_Q_skew_nz)
                for (index, (i, q)) in enumerate(row)
                    if i > j
                        θ = cube_space[i, j]
                        nij = q
                        FH_ij = norm(nij) * FH_ij_storage[i, j]
                        FL_ij = norm(nij) * FL_ij_storage[i, j]
                        du[i] += FH_ij + θ * FL_ij
                        du[j] -= FH_ij + θ * FL_ij
                    end
                end
            end
        else
            for (j, row) in enumerate(weird_Q_skew_nz)
                for (index, (i, q)) in enumerate(row)
                    if i > j
                        θ = cube_space[i, j]
                        nij = q
                        FH_ij = norm(nij) * FH_ij_storage[i, j]
                        FL_ij = norm(nij) * FL_ij_storage[i, j]
                        du[i] += θ * FL_ij + (1 - θ) * FH_ij
                        du[j] -= θ * FL_ij + (1 - θ) * FH_ij
                    end
                end
            end
        end
    end

    if entropy_inequality == :none || entropy_inequality == :local
        if blending_strat == :fft && filter_strength > 0
            store_fourier_coeffs!(Rdr, r_L .- r_H)
            lower_bounds = 1 .- exp.(-filter_strength * (1:length(Rdr)))
            Rdr .*= lower_bounds
            add_inverse_fourier!(r_H, Rdr)
        end
        
        du .= -(M \ r_H)
    elseif entropy_inequality == :semi_local
        # The only difference here is that r_H is stored in du instead
        if blending_strat == :fft && filter_strength > 0
            store_fourier_coeffs!(Rdr, r_L .- du)
            lower_bounds = 1 .- exp.(-filter_strength * (1:length(Rdr)))
            Rdr .*= lower_bounds
            add_inverse_fourier!(du, Rdr)
        end
        du .= -(M \ du)
    elseif entropy_inequality == :global && blend == :scalar
        # Simply take a convex combination (scalar lambda) between r_L and r_H to get an entropy stable method
        λ = clamp(v'r_L - (psi(u[end]) - psi(u[1])) / (v' * (r_L - r_H)), 0, 1)

        du .= -(M \ (λ * r_H + (1 - λ) * r_L))
    elseif entropy_inequality == :global && blend == :knapsack
        # Set blending coefficients as knapsack solution which satisfies a global entropy inequality
        # First, efficiently compute Δ'v (store in Dv)
        if blending_strat == :nodal
            Dv[1] = -v[1]
            for i in 1:(length(v) - 1)
                Dv[i + 1] = v[i] - v[i + 1]
            end
            Dv[end] = v[end]

            # Efficiently compute R (r^L - r^H) (store in Rdr)
            cum_sum = zero(eltype(Rdr))
            for i in eachindex(du)
                Rdr[i] = cum_sum
                cum_sum += r_L[i] - r_H[i]
            end
            Rdr[end] = zero(eltype(Rdr))
        elseif blending_strat == :fft
            # Handling SVectors is odd
            store_fourier_coeffs!(Dv, v)
            store_fourier_coeffs!(Rdr, r_L .- r_H)
        else
            throw("Blending strat $blending_strat not supported")
        end

        # Determine lower bounds
        lower_bounds = 1 .- exp.(-filter_strength * (1:length(a)))

        # Set a as the elementwise product between R (r^L - r^H) (in Rdr) and Δ'v (in Dv)
        a .= dot.(Rdr, Dv)
        b = -v'r_H - (psi(u[end]) - psi(u[1])) - a'lower_bounds

        # Now, run the knapsack solver
        knapsack_solver!(θ, a, b, upper_bounds=1 .- lower_bounds)

        # Finally, set the result
        Rdr .*= θ .+ lower_bounds
        if blending_strat == :nodal
            for i in eachindex(r_H)
                r_H[i] += Rdr[i + 1] - Rdr[i]
            end
        elseif blending_strat == :fft
            add_inverse_fourier!(r_H, Rdr)
        else
            throw("Blending strat $blending_strat not supported")
        end

        # Now, r_H stores r
        du .= -(M \ r_H)
    else
        throw("Not supported: blend: $blend, entropy_inequality: $entropy_inequality")
    end

    if !isnothing(bc)
        du[1] = zero(eltype(du))
        du[end] = zero(eltype(du))
    end

    return du
end

if is_periodic
    D = periodic_derivative_operator(; derivative_order=1, accuracy_order=accuracy_order, xmin=xmin, xmax=xmax, N=num_nodes)
else
    D = derivative_operator(WilliamsDuru2024(:central), derivative_order=1, accuracy_order=accuracy_order,
                               xmin=xmin, xmax=xmax, N=num_nodes) # works well: WilliamsDuru2024(:central)
end

M = mass_matrix(D)
Q = M * Matrix(D)
Q_skew = Q - Q'

x = SummationByPartsOperators.grid(D)

# Determine internal blending length
blending_length = length(x)

if blending_strat == :nodal
    blending_length = length(x) + 1
end

knapsack_solver = knapsack(blending_length)

# Construct an index pattern for Q_skew...
Q_skew_nz = zip(findnz(sparse(Q_skew))...)
weird_Q_skew_nz = Vector{Vector{Tuple{Int64, Float64}}}(undef, size(Q_skew, 2))

for j in 1:size(Q_skew, 2)
    weird_Q_skew_nz[j] = Vector{Tuple{Int64, Float64}}[]
end

for (i, j, q) in Q_skew_nz
    push!(weird_Q_skew_nz[j], (i, q))
end

# Which allows looping like this!
# Q_skew_copy = deepcopy(Q_skew)
# for (j, row) in enumerate(weird_Q_skew_nz)
#     for (index, (i, q)) in enumerate(row)
#         @assert Q_skew[i, j] == q
#         @assert Q_skew[j, i] == -q

#         Q_skew_copy[i, j] = 0.
#         Q_skew_copy[j, i] = 0.
#     end
# end

# @assert iszero(Q_skew_copy)