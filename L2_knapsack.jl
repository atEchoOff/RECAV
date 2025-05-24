struct QuadraticKnapsackMinimizer{Ttol}
    tol::Ttol
    a_over_w::Vector{Float64}

    function QuadraticKnapsackMinimizer{Ttol}(size::Int64; tol::Ttol = 100 * eps()) where {Ttol}
        return new{Ttol}(tol, zeros(Ttol, size))
    end
end

function (s::QuadraticKnapsackMinimizer)(x, a, b; upper_bounds=1., w=1., tol = 100 * eps(), maxit=200)
    # maxit is a huge upper bound here. The iteration will take at most N + 1 iterations, but usually takes around 1 - 3, and rarely 4.
    # Note also, if maxit is reached, it likely implies that b is negative, so the problem needs cleaning

    # Here we minimize theta subject to a'theta >= b

    if b <= 0.0
        # x = 0 is the optimal solution for the minimization
        # meaning the optimal solution for the maximization problem is upper_bounds
        fill!(x, zero(eltype(x)))
        return x
    end

    # Now, we want a'theta = b > 0. 

    # Suppose we find theta. Then, (-a)'theta = (-b), or a'theta = b. Then we are done! Proceed as normal

    # Start the Newton iteration
    lambdak = zero(eltype(a))
    itercount = 0 # for sanity check

    s.a_over_w .= a ./ w
    
    for _ in range(0, maxit)
        # Clip current solution within feasible domain
        x .= lambdak * s.a_over_w
        x .= clamp.(x, 0., upper_bounds) # FIXME slow

        f_val = a' * x - b

        if abs(f_val) / max(1, norm(a)) < tol
            break
        end

        deriv = a' * (s.a_over_w .* (x .< upper_bounds) .* (lambdak * s.a_over_w .>= 0.0) .* (s.a_over_w .> 0.0))

        # Compute the next root
        lambdak -= f_val / deriv

        itercount += 1 # for sanity check
    end

    ### These are all my sanity checks. Non well-posed problems may break them, so if issues are found, uncomment these and the sanity check comments above for checking.
    if itercount > 6
        println("The itercount was $itercount")
        @show a
        @show b
        @show upper_bounds
        println(a'x - b)
    end

    return x
end