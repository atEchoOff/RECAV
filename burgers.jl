xmin = -1
xmax = 1
is_periodic = true
T = 10.

include("common.jl")

equations = InviscidBurgersEquation1D()
initial_condition = initial_condition_burgers_gaussian

u0 = initial_condition.(x)
# @. u0 = prim2cons.(u0, equations)

include("initialize_globals.jl")

psi(u) = (1/6 * u .^ 3)[1]

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
    bc = [u0[1], u0[end]],
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
    plot(x, getindex.(sol.u[i], 1), leg=false, ylims=(0., 1.))
end