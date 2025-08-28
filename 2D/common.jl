using Polyester: @batch

function rhs!(du, u, cache, t)
    (; Q_skew, Q_skew_rows, Q_skew_vals, M, psi, alpha, dt, blend, entropy_inequality, low_order_volume_flux, equations, r_H, a, θ, v, knapsack_solvers, FH_ij_storage, FL_ij_storage, index_of_ji, flux_storage_1, flux_storage_2, flux) = cache

    @. v = cons2entropy.(u, equations)
    fill!(r_H, zero(eltype(r_H)))
    fill!(du, zero(eltype(du)))

    @Threads.threads :static for j in axes(Q_skew, 2)
        # Precompute fluxes
        flux_storage_1[j] = flux(u[j], 1, equations)
        flux_storage_2[j] = flux(u[j], 2, equations)
    end

    # This can be threaded since everything inside the loop should be completely indepedent
    @Threads.threads :static for j in axes(Q_skew, 2)
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
            end
        end

        # Part 2
        # Now, we compute the fluxes, and a, b (if necessary)
        for (index, skew_internal) in enumerate(col)
            i = Q_skew_rows[skew_internal]
            nij = Q_skew_vals[skew_internal]

            FH_ij_local[index] = FH_ij = (flux_storage_1[i] + flux_storage_1[j]) / 2 * nij[1] + (flux_storage_2[i] + flux_storage_2[j]) / 2 * nij[2]

            if entropy_inequality == :semi_local
                Ψij = (psi(u[j], nij) - psi(u[i], nij))

                if blend == :knapsack
                    FL_ij_local[index] = FL_ij = norm(nij) * low_order_volume_flux(u[i], u[j], nij / norm(nij), equations)
                    a_local[index] = dot(v[i] - v[j], FL_ij - FH_ij)
                elseif blend == :viscosity
                    FL_ij_local[index] = FL_ij = norm(nij) * (u[i] - u[j])
                    a_local[index] = dot(v[i] - v[j], FL_ij)
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
            if blend == :knapsack
                knapsack_solver_local!(θ_local, a_local, b_local)
            elseif blend == :viscosity
                squared_norm = dot(a_local, a_local)
                if b_local < 100 * eps() || squared_norm < 100 * eps()
                    # HOM already satisfies semi_local entropy
                    θ_local .= zero(eltype(θ_local))
                else
                    θ_local .= a_local
                    θ_local .*= b_local / squared_norm
                end
            end
        end
    end

    # Part 4
    # Blend together!
    if entropy_inequality == :semi_local
        for j in axes(Q_skew, 2)
            col = nzrange(Q_skew, j)

            FH_ij_local = @view FH_ij_storage[1:length(col), j]
            FL_ij_local = @view FL_ij_storage[1:length(col), j]
            θ_local = @view θ[1:length(col), j]

            for (index, skew_internal) in enumerate(col)
                i = Q_skew_rows[skew_internal]
                if i > j
                    FH_ij = FH_ij_local[index]
                    FL_ij = FL_ij_local[index]

                    θij = θ_local[index]
                    θji = θ[index_of_ji[i, j], i]

                    θij = max(θij, θji)

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

    return du
end

if is_periodic
    D = periodic_derivative_operator(; derivative_order=1, accuracy_order=accuracy_order, xmin=xmin, xmax=xmax, N=num_nodes)
else
    D = derivative_operator(WilliamsDuru2024(:central), derivative_order=1, accuracy_order=accuracy_order,
                               xmin=xmin, xmax=xmax, N=num_nodes) # works well: WilliamsDuru2024(:central)
end

M1d = mass_matrix(D)(num_nodes)
Q1d = sparse(M1d * Matrix(D))

M = kron(M1d, M1d)
Qx = kron(M1d, Q1d)
Qy = kron(Q1d, M1d)

Q = SVector{2, Float64}.(Qx, Qy)

Qt = permutedims(Q)

Q_skew = Q - Qt

x1d = SummationByPartsOperators.grid(D)

x = vec(@. x1d + 0 * x1d')
y = vec(@. 0 * x1d + x1d')

X = SVector{2, Float64}.(x, y)

# Construct an index pattern for Q_skew...
Q_skew_nz = zip(findnz(sparse(Q_skew))...)

weird_Q_skew_nz = Vector{Vector{Tuple{Int64, SVector{2, Float64}}}}(undef, size(Q_skew, 2))

for j in 1:size(Q_skew, 2)
    weird_Q_skew_nz[j] = Vector{Tuple{Int64, SVector{2, Float64}}}[]
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
    push!(knapsack_solvers, knapsack(length(col)))
end

# Finally,
Q_skew = sparse(Q_skew)
Q_skew_rows = rowvals(Q_skew)
Q_skew_vals = nonzeros(Q_skew)

# Do the same for Q' to get first nonzero entry in each row
Q_skew_T = sparse(permutedims(Q_skew))
Q_skew_T_row_vals = rowvals(Q_skew_T)

index_of_ji = spzeros(Int64, size(Q_skew))


for j in axes(Q_skew, 2)
    col = nzrange(Q_skew, j)
    for (index, k) in enumerate(col)
        i = Q_skew_rows[k]
        q = Q_skew_vals[k]

        rows_in_col_i = Q_skew_rows[nzrange(Q_skew, i)]
        
        index_of_ji[i, j] = findfirst(isequal(j), rows_in_col_i)
    end
end

# Test
# for j in axes(Q_skew, 2)
#     col = nzrange(Q_skew, j)
#     for (index, k) in enumerate(col)
#         i = Q_skew_rows[k]
#         q = Q_skew_vals[k]

#         @assert Q_skew[i, j] == q

#         FH_ij_storage[1, j] = SVector{2, Float64}(j, j)
#     end
# end

# for j in axes(Q_skew, 2)
#     col = nzrange(Q_skew, j)
#     for (index, k) in enumerate(col)
#         i = Q_skew_rows[k]
#         q = Q_skew_vals[k]

#         @assert Q_skew[i, j] == q
#         @assert Q_skew[j, i] == -q
#         @assert FH_ij_storage[1, i] == SVector{2, Float64}(i, i)
#     end
# end
