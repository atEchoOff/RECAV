using MAT

xmin = -5.
xmax = 5.
is_periodic = false
reflective_bcs = false
T = 1.8

include("common.jl")

equations = CompressibleEulerEquations1D(1.4)
initial_condition = initial_condition_shu_osher

u0 = initial_condition.(x)

include("initialize_globals.jl")

psi(u, nij) = u[2] * nij

cache = (;
    M, 
    psi, 
    preserve_positivity,
    dt,
    blend,
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
    bc = [u0[1], u0[end]],
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

weno_sol = matread("1D/weno5_shuosher.mat")
# plot(rd.Vp * md.x, u_plot, leg=false)
plot(weno_sol["x"][1:5:end], weno_sol["rho"][1:5:end], label="WENO", w=2)
plot!(x, getindex.(sol.u[end], 1), label="KL-FD", w=2)
plot!(xlim=(-3, 2.5))
plot!(ylim=(2.5, 5))