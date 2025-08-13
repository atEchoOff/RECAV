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
    (; Q_skew, Q_skew_rows, Q_skew_vals, M, psi, preserve_positivity, dt, potential_blend, blend, entropy_inequality, low_order_volume_flux, equations, r_H, r_H_temp, r_L, a, θ, v, knapsack_solvers, bc, is_periodic, weak_bcs, reflective_bcs, FH_ij_storage, FL_ij_storage, index_of_ji, flux_storage, flux, l_c, b_global) = cache

    @. v = cons2entropy.(u, equations)
    fill!(r_H, zero(eltype(r_H)))
    fill!(du, zero(eltype(du)))

    @static if preserve_positivity >= 0.
        fill!(r_H_temp, zero(eltype(r_H_temp)))
        fill!(r_L, zero(eltype(r_L)))
    end

    @static if !is_periodic && reflective_bcs && weak_bcs
        # Setup reflective BCs
#         Compute rho, v1, p from uBC and cons2prim
#         Set bc = prim2cons(SVector(rho, -v1, p), equations)
        bc = zeros(eltype(u0), 2)

        rho, v1, p = cons2prim(u[1], equations)
        bc[1] = prim2cons(SVector{3, Float64}(rho, -v1, p), equations)

        rho, v1, p = cons2prim(u[end], equations)
        bc[2] = prim2cons(SVector{3, Float64}(rho, -v1, p), equations)
    end
    
    @static if !is_periodic && weak_bcs
        du[1] += 1 / M[1, 1] * low_order_volume_flux(bc[1], u[1], SVector{1, Float64}(1.), equations)
        du[end] += 1 / M[end, end] * low_order_volume_flux(bc[end], u[end], SVector{1, Float64}(-1.), equations)
    end

    for i in axes(Q_skew, 1)
        flux_storage[i] = flux(u[i], 1, equations)
    end

    # This can be threaded since everything inside the loop should be completely indepedent
    for j in axes(Q_skew, 2)
        col = nzrange(Q_skew, j)
        # Part 1
        # Allocate necessary variables
        # vertical (contiguous) memory slice for arrays
        FH_ij_local = @view FH_ij_storage[1:length(col), j]
        FL_ij_local = @view FL_ij_storage[1:length(col), j]

        b_local = 0. # FIXME explicit type

        @static if potential_blend == :free
            a_local = @view a[1:length(col) + 1, j]
            θ_local = @view θ[1:length(col) + 1, j]
            a_local[end] = zero(eltype(a_local))
        else
            a_local = @view a[1:length(col), j]
            θ_local = @view θ[1:length(col), j]
        end

        knapsack_solver_local! = knapsack_solvers[j]

        # Part 2
        # Now, we compute the fluxes, and a, b (if necessary)
        for (index, skew_internal) in enumerate(col)
            i = Q_skew_rows[skew_internal]
            nij = Q_skew_vals[skew_internal]

            FH_ij_local[index] = FH_ij = nij * (flux_storage[i] + flux_storage[j]) / 2

            @static if entropy_inequality == :semi_local
                Ψij = psi(u[j], nij) - psi(u[i], nij)

                @static if blend == :knapsack
                    FL_ij_local[index] = FL_ij = norm(nij) * low_order_volume_flux(u[i], u[j], SVector{1}(nij / norm(nij)), equations)
                    @static if potential_blend == :free
                        a_local[index] = dot(v[i], FL_ij - FH_ij)
                        a_local[end] -= dot(v[j], FL_ij - FH_ij)
                    else
                        a_local[index] = dot(v[i] - v[j], FL_ij - FH_ij)
                    end
                elseif blend == :viscosity
                    FL_ij_local[index] = FL_ij = norm(nij) * (u[i] - u[j])
                    @static if potential_blend == :free
                        a_local[index] = dot(v[i], FL_ij)
                        a_local[end] -= dot(v[j], FL_ij)
                    else
                        a_local[index] = dot(v[i] - v[j], FL_ij)
                    end
                end

                b_local += dot(v[j] - v[i], FH_ij) - Ψij

                @static if preserve_positivity >= 0.
                    # Build up r_L and r_H for positivity
                    if i > j
                        r_H_temp[i] += FH_ij
                        r_H_temp[j] -= FH_ij

                        r_L[i] += FL_ij
                        r_L[j] -= FL_ij
                    end
                end
            elseif blend == :LOM
                FL_ij_local[index] = norm(nij) * low_order_volume_flux(u[i], u[j], SVector{1}(nij / norm(nij)), equations)
            end
        end

        # Part 3
        # We determine the blending coefficients (if necessary)
        @static if entropy_inequality == :semi_local
            @static if blend == :knapsack && preserve_positivity < 0
                # If we are preserving positivity, we are not ready to call this yet. We need upper bounds
                knapsack_solver_local!(θ_local, a_local, b_local)
            elseif blend == :knapsack && preserve_positivity >= 0
                # We should save b_local into b_global
                b_global[j] = b_local
            elseif blend == :viscosity
                norm_squared = dot(a_local, a_local)
                if b_local < 100 * eps() || norm_squared < 100 * eps()
                    # HOM already satisfies semi_local entropy
                    θ_local .= zero(eltype(θ_local))
                else
                    θ_local .= a_local
                    θ_local .*= b_local / norm_squared
                end
            end
        end
    end

    @static if entropy_inequality == :semi_local && preserve_positivity >= 0
        # Compute positivity constants first
        # Set positions in r_H_temp
        for i in axes(Q_skew, 1)
            bound_1 = 0.
            if r_H_temp[i][1] - r_L[i][1] > 100 * eps()
                bound_1 = 1 - (1 - preserve_positivity) / dt * (M[i, i] * u[i][1] - dt * r_L[i][1]) / (r_H_temp[i][1] - r_L[i][1])
            end

            bound_3 = 0.
            if r_H_temp[i][3] - r_L[i][3] > 100 * eps()
                bound_3 = 1 - (1 - preserve_positivity) / dt * (M[i, i] * u[i][3] - dt * r_L[i][3]) / (r_H_temp[i][3] - r_L[i][3])
            end

            r_H_temp[i] = SVector{3, Float64}(bound_1, 0., bound_3)
        end

        # The blend better be knapsack or this isnt going to work
        @static if blend !== :knapsack
            println("WARNING: for positivity preservation, blend is overridden to knapsack")
        end

        for j in axes(Q_skew, 2)
            col = nzrange(Q_skew, j)

            @static if potential_blend == :free
                θ_local = @view θ[1:length(col) + 1, j]
                l_c_local = @view l_c[1:length(col) + 1, j]
                a_local = @view a[1:length(col) + 1, j]

                l_c_local[end] = 0.
            else
                θ_local = @view θ[1:length(col), j]
                l_c_local = @view l_c[1:length(col), j]
                a_local = @view a[1:length(col), j]
            end

            b_local = b_global[j]
            knapsack_solver_local! = knapsack_solvers[j]

            for (index, skew_internal) in enumerate(col)
                i = Q_skew_rows[skew_internal]
                l_c_local[index] = max(r_H_temp[i][1], r_H_temp[i][3], r_H_temp[j][1], r_H_temp[j][3], 0)
            end

            # Adjust constants for positivity
            b_local -= dot(a_local, l_c_local)

            knapsack_solver_local!(θ_local, a_local, b_local; upper_bounds=1 .- l_c_local)
            θ_local .+= l_c_local
        end
    end

    for j in axes(Q_skew, 2)
        col = nzrange(Q_skew, j)

        FH_ij_local = @view FH_ij_storage[1:length(col), j]
        FL_ij_local = @view FL_ij_storage[1:length(col), j]
        θ_local = @view θ[1:length(col), j]

        # Part 4
        # Blend together!
        for (index, skew_internal) in enumerate(col)
            i = Q_skew_rows[skew_internal]

            if i > j
                FH_ij = FH_ij_local[index]
                FL_ij = FL_ij_local[index]

                θij = θ_local[index]
                θji = θ[index_of_ji[i, j], i]
                θij = max(θij, θji)

                # The rest is simple!
                @static if blend == :viscosity
                    r_H[i] += FH_ij + θij * FL_ij
                    r_H[j] -= FH_ij + θij * FL_ij
                elseif blend == :knapsack
                    r_H[i] += θij * FL_ij + (1 - θij) * FH_ij
                    r_H[j] -= θij * FL_ij + (1 - θij) * FH_ij
                elseif blend == :HOM
                    r_H[i] += FH_ij
                    r_H[j] -= FH_ij
                elseif blend == :LOM
                    r_H[i] += FL_ij
                    r_H[j] -= FL_ij
                end
            end
        end
    end

    for i in eachindex(du)
        du[i] += -r_H[i] / M[i, i]
    end

    # Finish off boundary conditions
    if !is_periodic && !weak_bcs
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

max_length = 0
knapsack_solvers = knapsack[]
for (j, col) in enumerate(weird_Q_skew_nz)
    global max_length # wtf lol
    if length(col) > max_length
        max_length = length(col)
    end
    push!(knapsack_solvers, knapsack(length(col) + 1))
end

# Finally,
Q_skew = sparse(Q_skew)
Q_skew_rows = rowvals(Q_skew)
Q_skew_vals = nonzeros(Q_skew)

index_of_ji = spzeros(Int64, size(Q_skew))

for j in axes(Q_skew, 2)
    col = nzrange(Q_skew, j)
    for (index, k) in enumerate(col)
        local i
        i = Q_skew_rows[k]
        q = Q_skew_vals[k]

        rows_in_col_i = Q_skew_rows[nzrange(Q_skew, i)]
        
        index_of_ji[i, j] = findfirst(isequal(j), rows_in_col_i)
    end
end