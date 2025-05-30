struct QuadraticKnapsackMaximizer{Ttol}
    tol::Ttol
    a_over_w::Vector{Float64}
    direction::Int64

    function QuadraticKnapsackMaximizer{Ttol}(size::Int64;tol::Ttol = 100 * eps()) where {Ttol}
        return new{Ttol}(tol, zeros(Ttol, size), -1)
    end
end

function (s::QuadraticKnapsackMaximizer)(x, a, b; upper_bounds=1., w=1., tol = 100 * eps(), maxit=200)
    # maxit is a huge upper bound here. The iteration will take at most N + 1 iterations, but usually takes around 1 - 3, and rarely 4.
    # Note also, if maxit is reached, it likely implies that b is negative, so the problem needs cleaning

    # Use local knapsack memory
    a_over_w = @view s.a_over_w[1:length(a)]

    # Problem infeasibility... this is bad
    # Just return low order solution... best we can do
    if b <= 0.0
        fill!(x, zero(eltype(x)))
        # println("Infeasibility Detected")
        return x
    end

    # original_b = b # for sanity check
    b = sum(a * upper_bounds) - b

    if b <= 0.0
        # x = 0 is the optimal solution for the minimization
        # meaning the optimal solution for the maximization problem is upper_bounds
        x .= upper_bounds
        return x
    end

    # Start the Newton iteration
    lambdak = zero(eltype(a))
    itercount = 0 # for sanity check

    a_over_w .= a ./ w
    
    for _ in range(0, maxit)
        # Clip current solution within feasible domain
        x .= lambdak * a_over_w
        x .= clamp.(x, 0.0, upper_bounds)

        f_val = a' * x - b

        if abs(f_val) / max(1, norm(a)) < tol
            break
        end

        deriv = a' * (a_over_w .* (x .< upper_bounds) .* (lambdak * a_over_w .>= 0.0) .* (a_over_w .> 0.0))

        # Compute the next root
        lambdak -= f_val / deriv

        itercount += 1 # for sanity check
    end

    # x .= sgn_a .* x

    x .= upper_bounds .- x

    ### These are all my sanity checks. Non well-posed problems may break them, so if issues are found, uncomment these and the sanity check comments above for checking.
    # if itercount > 3
    #     println("The itercount was $itercount")
    #     @show a
    #     @show b
    #     @show upper_bounds
    #     @show x
    #     println(a'x - b)

    #     sleep(1)
    # end

    # @assert all(x .== 1.0)

    # @assert itercount <= 4

    # @assert a' * x <= original_b + tol
    # @assert abs(a' * x - original_b) <= tol
    # @assert all(x .<= upper_bounds)
    # @assert all(x .>= 0.0)
    return x
end