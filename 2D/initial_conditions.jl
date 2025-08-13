function initial_condition_density_wave_fast(x, y, t)
    u = 1.7
    v = 1.
    rho = .5 * sin(pi * (x - 1.7t)) * sin(pi * (y - t)) + 1
    p = 1.

    return SVector(prim2cons(SVector(rho, u, v, p), equations))
end

function initial_condition_football(x, y, t)
    if (0 <= x) & (0 <= y)
        rho = .5313
        u, v = 0, 0
        p = 0.4
    elseif (x < 0) & (0 <= y)
        rho = 1
        u, v = .7276, 0
        p = 1
    elseif (x < 0) & (y < 0)
        rho = .8
        u, v = 0, 0
        p = 1
    elseif (0 <= x) & (y < 0)        
        rho = 1
        u, v = 0, .7276
        p = 1
    else 
        @show x, y
    end

    return prim2cons(SVector(rho, u, v, p), equations)
end

function initial_condition_khi(x, y, t)
    # KHI
    B = tanh(15 * y + 7.5) - tanh(15 * y - 7.5)
    rho = 0.5 + 0.75 * B
    u = 0.5 * (B - 1)
    v = 0.1 * sin(2 * pi * x)
    p = 1
    return prim2cons(SVector(rho, u, v, p), equations)
end

function initial_condition_khi_competitive(x, y, t)
    # KHI
    if (y > .25)
        rho = 1
        u = .5
        v = .01sin(2pi * x)
        p = 2.5
    elseif (y < -.25)
        rho = 1
        u = .5
        v = .01sin(2pi * x)
        p = 2.5
    else
        rho = 2
        u = -.5
        v = .01sin(2pi * x)
        p = 2.5
    end

    return prim2cons(SVector(rho, u, v, p), equations)
end

function initial_condition_BM_vortex(x, y, t)
    pbar = 9.0 / equations.gamma
    delta = 0.05
    epsilon = 30
    H = (y < 0) ? tanh(epsilon * (y + 0.25)) : tanh(epsilon * (0.25 - y))
    rho = 1.0
    v1 = H
    v2 = delta * cos(2.0 * pi * x)
    p = pbar
    return prim2cons(SVector(rho, v1, v2, p), equations)
end

function initial_condition_shear_rollup(x, y, t)
    p = 9.0 / equations.gamma
    rho = pi / 15
    delta = .05

    if y <= pi
        u = tanh((y - pi/2)/rho)
    else
        u = tanh((3pi/2 - y)/rho)
    end

    v = delta * sin(x)

    return prim2cons(SVector(rho, u, v, p), equations)
end