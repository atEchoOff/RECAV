xmin = -1
xmax = 1
is_periodic = true
T = 2.

include("common.jl")

equations = LinearScalarAdvectionEquation1D(1.)
initial_condition = initial_condition_advection_sin

u0 = initial_condition.(x, 0.)
# @. u0 = prim2cons.(u0, equations)

include("initialize_globals.jl")

psi(u) = .5 * u^2

cache = (;
    M, 
    psi, 
    alpha = preserve_positivity,
    dt,
    blend,
    entropy_inequality, 
    volume_flux, 
    low_order_volume_flux,
    equations, 
    r_H, 
    a, 
    θ, 
    v,
    knapsack_solvers,
    bc = nothing,
    weak_bcs,
    Q_skew,
    Q_skew_rows,
    Q_skew_vals,
    FH_ij_storage,
    FL_ij_storage
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