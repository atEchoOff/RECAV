xmin = -1
xmax = 1
is_periodic = false
reflective_bcs = false
T = 1.

include("common.jl")

equations = InviscidBurgersEquation1D()
initial_condition = initial_condition_burgers_gaussian

u0 = initial_condition.(x)
# @. u0 = prim2cons.(u0, equations)

include("initialize_globals.jl")

psi(u, nij) = (1/6 * u .^ 3)[1] * nij

cache = (;
    M, 
    psi, 
    preserve_positivity,
    dt,
    blend,
    potential_blend,
    entropy_inequality, 
    volume_flux, 
    low_order_volume_flux,
    equations, 
    r_H, 
    r_H_temp,
    r_L,
    a, 
    θ, 
    v,
    knapsack_solvers,
    bc = nothing,
    is_periodic,
    weak_bcs,
    reflective_bcs,
    Q_skew,
    Q_skew_rows,
    Q_skew_vals,
    FH_ij_storage,
    FL_ij_storage,
    index_of_ji, 
    flux_storage, 
    flux,
    l_c,
    b_global
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