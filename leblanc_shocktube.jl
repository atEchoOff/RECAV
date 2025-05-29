xmin = -10
xmax = 10
is_periodic = false
T = 1e-4

include("common.jl")

equations = CompressibleEulerEquations1D(1.4)
initial_condition = initial_condition_leblanc_shocktube

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
            callback=AliveCallback(alive_interval=100), 
            adaptive=adaptive)

u = cons2prim.(sol.u[end], equations)

plot(x, getindex.(u, 2), lw=2)