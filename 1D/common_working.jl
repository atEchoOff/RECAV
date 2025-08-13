import Trixi: flux_hllc
# use the 2D implementation since the 1D version doesn't account for n::Integer < 0
function flux_hllc(u_ll, u_rr, n::SVector{1}, equations::CompressibleEulerEquations1D)
    f = flux_hllc(SVector(u_ll[1], u_ll[2], 0, u_ll[3]), 
                  SVector(u_rr[1], u_rr[2], 0, u_rr[3]), 
                  SVector(n[1], 0.0), 
                  CompressibleEulerEquations2D(equations.gamma))    
    return SVector(f[1], f[2], f[4])
end

function rhs!(du, u, cache, t)
    (; Q_skew, Q_skew_rows, Q_skew_vals, M, psi, alpha, dt, blend, entropy_inequality, volume_flux, low_order_volume_flux, equations, r_H, a, θ, v, knapsack_solvers, bc, weak_bcs, FH_ij_storage, FL_ij_storage) = cache

    @. v = cons2entropy.(u, equations)
    fill!(r_H, zero(eltype(r_H)))
    fill!(du, zero(eltype(du)))

    if !isnothing(bc) && weak_bcs
        du[1] += 1 / M[1, 1] * low_order_volume_flux(bc[1], u[1], SVector{1, Float64}(1.), equations)
        du[end] += 1 / M[end, end] * low_order_volume_flux(bc[end], u[end], SVector{1, Float64}(-1.), equations)
    end

    # This can be threaded since everything inside the loop should be completely indepedent
    for j in axes(Q_skew, 2)
        col = nzrange(Q_skew, j)
        # Part 1
        # Allocate necessary variables
        # vertical (contiguous) memory slice for arrays
        FH_ij_local = @view FH_ij_storage[1:length(col), j]
        FL_ij_local = @view FL_ij_storage[1:length(col), j]

        if entropy_inequality == :semi_local
            b_local = 0. # FIXME explicit type

            a_local = @view a[1:length(col), j]
            θ_local = @view θ[1:length(col), j]
            if blend == :knapsack
                knapsack_solver_local! = knapsack_solvers[j]
            elseif blend == :viscosity || blend == :scalar
                # For scalar, a_local = sum(a). For visc, a_local = norm(a)^2. Either way, we sum. So, initialize to 0. 
                a_local[1] = zero(eltype(a_local))
            else
                println("Blending strat $blend not supported")
                return du
            end
        end

        # Part 2
        # Now, we compute the fluxes, and a, b (if necessary)
        for (index, skew_internal) in enumerate(col)
            i = Q_skew_rows[skew_internal]
            nij = Q_skew_vals[skew_internal]

            FH_ij_local[index] = FH_ij = abs(nij) * volume_flux(u[i], u[j], SVector{1}(nij / norm(nij)), equations)

            if entropy_inequality == :semi_local
                Ψij = (psi(u[j]) - psi(u[i])) * nij

                if blend == :knapsack
                    FL_ij_local[index] = FL_ij = abs(nij) * low_order_volume_flux(u[i], u[j], SVector{1}(nij / norm(nij)), equations)
                    a_local[index] = dot(v[i] - v[j], FL_ij - FH_ij)
                elseif blend == :scalar
                    FL_ij_local[index] = FL_ij = abs(nij) * low_order_volume_flux(u[i], u[j], SVector{1}(nij / norm(nij)), equations)
                    a_local[1] += dot(v[i] - v[j], FL_ij - FH_ij)
                elseif blend == :viscosity
                    FL_ij = abs(nij) * (u[i] - u[j])
                    # We use the analytical solution for viscosity for θ. Instead of being scalar, everyhing is multiplied by elements of a. So, we will take this into account
                    a_component = dot(v[i] - v[j], FL_ij)

                    FL_ij_local[index] = FL_ij = a_component * FL_ij
                    a_local[1] += a_component^2
                end

                b_local += dot(v[j] - v[i], FH_ij) - Ψij
            elseif entropy_inequality == :none
                # We are done with just the high order method
                if i > j
                    r_H[i] += FH_ij
                    r_H[j] -= FH_ij
                end
            else
                println("Entropy inequality $entropy_inequality is not supported")
                return du
            end
        end

        # Part 3
        # We determine the blending coefficients (if necessary)
        if entropy_inequality == :semi_local
            # TODO limiting coefficients. The computation of r_H and r_L can be done in previous loop, making this a bit easier
            if blend == :knapsack
                knapsack_solver_local!(θ_local, a_local, b_local)
            elseif blend == :scalar
                if b_local < 100 * eps() || abs(a_local[1]) < 100 * eps()
                    # HOM already satisfies semi_local entropy
                    θ_local[1] = zero(eltype(θ_local))
                else
                    θ_local[1] = clamp(b_local / a_local[1], 0., 1.)
                end
            elseif blend == :viscosity
                # This is a weird one. FL_ij is really a_i * FL_ij. So, rather than θ = ba / a'a, its θ = b / a'a = b / a_local[1]
                if b_local < 100 * eps() || abs(a_local[1]) < 100 * eps()
                    # HOM already satisfies semi_local entropy
                    θ_local[1] = zero(eltype(θ_local))
                else
                    θ_local[1] = b_local / a_local[1] # no clamping necessary, a_local > 100 * eps() for SURE for sure
                    a[2] += b_local
                    a[3] += 1
                end
            end
        end

        # Part 4
        # Blend together!
        if entropy_inequality == :semi_local
            for (index, skew_internal) in enumerate(col)
                i = Q_skew_rows[skew_internal]
                if i > j
                    FH_ij = FH_ij_local[index]
                    FL_ij = FL_ij_local[index]

                    # How we get blending coeffs depends on blend
                    if blend == :knapsack
                        θij = θ_local[index]
                    else
                        θij = θ_local[1]
                    end

                    # The rest is simple!
                    if blend == :viscosity
                        r_H[i] += FH_ij + θij * FL_ij
                        r_H[j] -= FH_ij + θij * FL_ij
                    else
                        r_H[i] += θij * FL_ij + (1 - θij) * FH_ij
                        r_H[j] -= θij * FL_ij + (1 - θij) * FH_ij
                    end
                end
            end
        end
    end

    for i in eachindex(du)
        du[i] += -r_H[i] / M[i, i]
    end

    # Finish off boundary conditions
    if !isnothing(bc) && !weak_bcs
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

# Construct an index pattern for Q_skew...
Q_skew_nz = zip(findnz(sparse(Q_skew))...)
weird_Q_skew_nz = Vector{Vector{Tuple{Int64, Float64}}}(undef, size(Q_skew, 2))

for j in 1:size(Q_skew, 2)
    weird_Q_skew_nz[j] = Vector{Tuple{Int64, Float64}}[]
end

for (i, j, q) in Q_skew_nz
    push!(weird_Q_skew_nz[j], (i, q))
end

knapsack_solvers = knapsack[]
for (j, col) in enumerate(weird_Q_skew_nz)
    push!(knapsack_solvers, knapsack(length(col)))
end

# Finally,
Q_skew = sparse(Q_skew)
Q_skew_rows = rowvals(Q_skew)
Q_skew_vals = nonzeros(Q_skew)