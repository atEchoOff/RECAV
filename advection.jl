xmin = -1
xmax = 1
T = 1.

include("common.jl")

equations = LinearScalarAdvectionEquation1D(1.)
initial_condition = initial_condition_advection_buzz

u0 = initial_condition.(x, 0.)
# @. u0 = prim2cons.(u0, equations)

include("initialize_globals.jl")

psi(u) = .5 * u^2

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

@gif for i in eachindex(sol.t)
    plot(x, sol.u[i], leg=false, ylims=(0, 2))
    plot!(x, initial_condition.(x, sol.t[i]))
end