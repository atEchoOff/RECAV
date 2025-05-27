using MAT

xmin = -5
xmax = 5
T = 1.8

include("common.jl")

equations = CompressibleEulerEquations1D(1.4)
initial_condition = initial_condition_shu_osher

u0 = initial_condition.(x)

include("initialize_globals.jl")

psi(u) = u[2]

cache = (;
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
    knapsack_solver! = knapsack_solver,
    bc = [u0[1] u0[end]],
    Q_skew_nz = Q_skew_nz,
    weird_Q_skew_nz
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

weno_sol = matread("weno5_shuosher.mat")
# plot(rd.Vp * md.x, u_plot, leg=false)
plot(weno_sol["x"][1:5:end], weno_sol["rho"][1:5:end], label="WENO", w=2)
plot!(x, getindex.(sol.u[end], 1), label="DG", w=2)