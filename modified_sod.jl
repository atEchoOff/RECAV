xmin = 0
xmax = 1
is_periodic = false
T = .2

include("common.jl")

equations = CompressibleEulerEquations1D(1.4)
initial_condition = initial_condition_modified_sod

u0 = initial_condition.(x)

include("initialize_globals.jl")

psi(u) = u[2]

cache = (;
    M, 
    psi, 
    alpha = preserve_positivity,
    dt,
    blend,
    entropy_inequality, 
    entropy_blend,
    blending_strat,
    filter_strength,
    volume_flux, 
    low_order_volume_flux,
    equations, 
    r_H, 
    r_L, 
    r_entropy_rhs, 
    a, 
    θ, 
    v,
    Rdr,
    Dv,
    knapsack_solver! = knapsack_solver,
    bc = [u0[1] u0[end]],
    Q_skew_nz = Q_skew_nz,
    weird_Q_skew_nz,
    cube_space
)

ode = ODEProblem(rhs!, u0, (0., T), cache)
sol = solve(ode, 
            timestepper, 
            dt = dt, 
            abstol=abstol, 
            reltol=reltol, 
            saveat=saveat, 
            callback=AliveCallback(alive_interval=1000), 
            adaptive=adaptive)

@gif for i in eachindex(sol.t)
    plot(x, getindex.(sol.u[i], 1), leg=false, ylims=(0, 1))
end

plot(x, getindex.(sol.u[end], 1), lw=2)