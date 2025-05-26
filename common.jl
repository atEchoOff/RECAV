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
    (; Q_skew_nz, M, psi, blend, blending_strat, filter_strength, volume_flux, equations, r_H, r_L, r_entropy_rhs, a, θ, v, knapsack_solver!, bc) = cache
    fill!(r_H, zero(eltype(du)))
    fill!(r_L, zero(eltype(du)))
    fill!(r_entropy_rhs, zero(eltype(x)))

    if !isnothing(bc)
        u[1] = bc[:, 1]
        u[end] = bc[:, end]
    end

    @. v = cons2entropy.(u, equations)

    for (i, j, q) in Q_skew_nz
        if i > j
            if abs(q) > 1e-12
                nij = q
                FH_ij = norm(nij) * volume_flux(u[i], u[j], SVector{1}(nij / norm(nij)), equations)
                FL_ij = norm(nij) * flux_lax_friedrichs(u[i], u[j], SVector{1}(nij / norm(nij)), equations)
                entropy_ij = norm(nij) * (v[j]'flux_central(u[i], u[j], SVector{1}(nij / norm(nij)), equations) - (psi(u[j]) - psi(u[i])) * nij / norm(nij))

                if blend == :local_entropy_local_scalar

                    λ = 0.
                    if abs(v[i]' * (FH_ij - FL_ij)) > 100 * eps() && -v[i]'FH_ij + entropy_ij > 100 * eps()
                        λ = clamp((v[i]'FH_ij - entropy_ij) / (v[i]' * (FH_ij - FL_ij)), 0., 1.)
                    end

                    r_H[i] += λ * FL_ij + (1 - λ) * FH_ij
                    r_H[j] -= λ * FL_ij + (1 - λ) * FH_ij
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
    end

    if blend == :highorder || blend == :local_entropy_local_scalar
        if blending_strat == :fft && filter_strength > 0
            store_fourier_coeffs!(Rdr, r_L .- r_H)
            lower_bounds = 1 .- exp.(-filter_strength * (1:length(a)))
            Rdr .*= lower_bounds
            add_inverse_fourier!(r_H, Rdr)
        end
        
        du .= -(M \ r_H)
    elseif blend == :loworder
        du .= -(M \ r_L)
    elseif blend == :global_entropy_global_scalar
        # Simply take a convex combination (scalar lambda) between r_L and r_H to get an entropy stable method
        λ = clamp(v'r_L - (psi(u[end]) - psi(u[1])) / (v' * (r_L - r_H)), 0, 1)

        du .= -(M \ (λ * r_H + (1 - λ) * r_L))
    elseif blend == :global_entropy_knapsack
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
        throw("Blending strategy $blend is not supported")
    end

    if !isnothing(bc)
        du[1] = zero(eltype(du))
        du[end] = zero(eltype(du))
    end

    return du
end

# D = periodic_derivative_operator(; derivative_order=1, accuracy_order=accuracy_order, xmin=xmin, xmax=xmax, N=num_nodes)
D = derivative_operator(WilliamsDuru2024(:central), derivative_order=1, accuracy_order=accuracy_order,
                               xmin=xmin, xmax=xmax, N=num_nodes)
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