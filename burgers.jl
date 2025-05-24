xmin = -1
xmax = 1
T = 10.

include("common.jl")

equations = InviscidBurgersEquation1D()
initial_condition = initial_condition_burgers_gaussian

u0 = initial_condition.(x)
# @. u0 = prim2cons.(u0, equations)

include("initialize_globals.jl")

psi(u) = (1/6 * u .^ 3)[1]

cache = (;
    Q_skew, 
    M, 
    psi, 
    blend,
    blending_strat,
    filter_strength,
    volume_flux, 
    equations, 
    r_H, 
    r_L, 
    r_entropy_rhs, 
    a, 
    θ, 
    v,
    Rdr,
    Dv,
    knapsack_solver! = knapsack_solver
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
#     plot(x, getindex.(sol_lom.u[i], 1), leg=false)
#     plot!(x, getindex.(sol_local.u[i], 1), leg=false)
# end