function initial_condition_advection_sin(x, t)
    return .5sin(2pi * (x - t)) + 1
end

function initial_condition_advection_buzz(x, t)
    return 1. * ceil(Int64, .5sin(pi * (x - t)) + 1)
end

function initial_condition_burgers_gaussian(x)
    return SVector(exp(-10 * x^2))
end

function initial_condition_modified_sod(x)
    if x[1] < .3
        rho = 1.0
        v1 = .75
        p = 1.0
    else
        rho = .125
        v1 = 0.0
        p = .1
    end
    return prim2cons(SVector(rho, v1, p), equations)
end

function initial_condition_shu_osher(x)
    if x[1] < -4
        rho = 3.857143
        v1 = 2.629369
        p = 10.3333 
    else
        rho = 1 + .2 * sin(5 * x[1])
        v1 = 0.0
        p = 1.0
    end
    return prim2cons(SVector(rho, v1, p), equations)
end

function initial_condition_density_wave_fast(x, t)
    u = 1.7
    rho = .5 * sin(pi * (x - 1.7t)) + 1
    p = 1.

    return SVector(prim2cons(SVector(rho, u, p), equations))
end