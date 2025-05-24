function rhs!(du, u, cache, t)
    (; Q_skew, M, psi, blend, blending_strat, filter_strength, volume_flux, equations, r_H, r_L, r_entropy_rhs, a, θ, v, knapsack_solver!) = cache
    fill!(r_H, zero(eltype(du)))
    fill!(r_L, zero(eltype(du)))
    fill!(r_entropy_rhs, zero(eltype(x)))

    @. v = cons2entropy.(u, equations)

    for i in eachindex(u), j in eachindex(u)
        if i > j
            if abs(Q_skew[i, j]) > 1e-12
                nij = Q_skew[i, j]
                FH_ij = norm(nij) * volume_flux(u[i], u[j], SVector{1}(nij / norm(nij)), equations)
                FL_ij = norm(nij) * flux_lax_friedrichs(u[i], u[j], SVector{1}(nij / norm(nij)), equations)
                entropy_ij = norm(nij) * (v[j]'flux_lax_friedrichs(u[i], u[j], SVector{1}(nij / norm(nij)), equations) - (psi(u[j]) - psi(u[i])) * nij / norm(nij))

                if blend == :local_entropy_local_scalar
                    # Low order method is entropy stable
                    # @assert -v[i]'FL_ij + entropy_ij <= 0

                    λ = 0.
                    if v[i]' * (FH_ij - FL_ij) > 100 * eps()
                        λ = clamp((v[i]'FH_ij - entropy_ij) / (v[i]' * (FH_ij - FL_ij)), 0., 1.)
                    end

                    r_H[i] += λ * FL_ij + (1 - λ) * FH_ij
                    r_H[j] -= λ * FL_ij + (1 - λ) * FH_ij
                else
                    r_H[i] += FH_ij
                    r_H[j] -= FH_ij

                    r_L[i] += FL_ij
                    r_L[j] -= FL_ij

                    r_entropy_rhs[i] += entropy_ij
                    r_entropy_rhs[j] -= entropy_ij
                end
            end
        end
    end

    if blend == :highorder || blend == :local_entropy_local_scalar
        du .= -(M \ r_H)
    elseif blend == :loworder
        # The low order method should be globally entropy stable
        @assert v'r_L >= 0

        du .= -(M \ r_L)

        # The low order method should also be nodally entropy stable
        r_entropy_rhs .= M \ r_entropy_rhs
        for i in eachindex(du)
            @assert v[i]'du[i] + r_entropy_rhs[i] <= 0
        end
    elseif blend == :global_entropy_global_scalar
        # Simply take a convex combination (scalar lambda) between r_L and r_H to get an entropy stable method
        λ = clamp(v'r_L / (v' * (r_L - r_H)), 0, 1)

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
            if eltype(u0) == Float64
                Dv .= real.(fft(v))
                Rdr .= real.(fft(r_L .- r_H))
            else
                # This is an SVector. fft each component
                Dv .= SVector{length(eltype(u0)), Float64}.([real.(fft(getindex.(v, i))) for i in 1:length(eltype(u0))]...)
                Rdr .= SVector{length(eltype(u0)), Float64}.([real.(fft(getindex.(r_L .- r_H, i))) for i in 1:length(eltype(u0))]...)
            end
        else
            throw("Blending strat $blending_strat not supported")
        end

        # Determine lower bounds
        lower_bounds = 1 .- exp.(-filter_strength * (1:length(a)))

        # Set a as the elementwise product between R (r^L - r^H) (in Rdr) and Δ'v (in Dv)
        a .= dot.(Rdr, Dv)
        b = -v'r_H - a'lower_bounds

        # Some naiive shock capturing...
        # b *= 10

        # Now, run the knapsack solver
        knapsack_solver!(θ, a, b, upper_bounds=1 .- lower_bounds)

        # Finally, set the result
        Rdr .*= θ .+ lower_bounds
        if blending_strat == :nodal
            for i in eachindex(r_H)
                r_H[i] += Rdr[i + 1] - Rdr[i]
            end
        elseif blending_strat == :fft
            if eltype(u0) == Float64
                r_H .+= real.(ifft(Rdr))
            else
                # We are working with SVectors again
                r_H .+= SVector{length(eltype(u0)), Float64}.([real.(ifft(getindex.(Rdr, i))) for i in 1:length(eltype(u0))]...)
            end
        else
            throw("Blending strat $blending_strat not supported")
        end

        # Now, r_H stores r
        du .= -(M \ r_H)
    else
        throw("Blending strategy $blend is not supported")
    end

    return du
end

D = periodic_derivative_operator(; derivative_order=1, accuracy_order=accuracy_order, xmin=xmin, xmax=xmax, N=num_nodes)
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