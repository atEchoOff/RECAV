function initial_condition_advection_sin(x, t)
    return .5sin(2pi * (x - t)) + 1
end

function initial_condition_advection_buzz(x, t)
    return 1. * ceil(Int64, .5sin(pi * (x - t)) + 1)
end

function initial_condition_burgers_gaussian(x)
    return SVector(exp(-10 * x^2))
end