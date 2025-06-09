xmin = -1
xmax = 1
is_periodic = true
T = 1

include("common.jl")

equations = CompressibleEulerEquations1D(1.4)
initial_condition = initial_condition_density_wave_fast

u0 = initial_condition.(x, 0.)

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
    # knapsack_solver! = knapsack_solver,
    knapsack_solvers = knapsack_solvers,
    bc = nothing,
    Q_skew_nz = Q_skew_nz,
    weird_Q_skew_nz,
    cube_space = cube_space,
    FH_ij_storage,
    FL_ij_storage,
    knapsack_shock_capturing
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

# @gif for i in eachindex(sol.t)
#     plot(x, getindex.(sol.u[i], 1), leg=false, ylims=(.5, 1.5))
#     plot!(x, getindex.(initial_condition.(x, sol.t[i]), 1), leg=false, ylims=(.5, 1.5))
# end

L2_error = sqrt(sum(M * map(x -> sum(x.^2), initial_condition.(x, sol.t[end]) - sol.u[end])))