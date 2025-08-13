function compute_vorticity(v1, v2)
    Dx = (M \ Qx)
    Dy = (M \ Qy)
    dv1dy = Dy * v1
    dv2dx = Dx * v2

    return dv2dx - dv1dy
end

xmin = -.5
xmax = .5
is_periodic = true
T = 1.

include("common.jl")

equations = CompressibleEulerEquations2D(1.4)
initial_condition = initial_condition_BM_vortex

u0 = initial_condition.(x, y, 0.)

include("initialize_globals.jl")

function psi(u, nij)
    rho, rho_v1, rho_v2, _ = u
    return dot(SVector(rho_v1, rho_v2), nij)
end

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
    FL_ij_storage,
    index_of_ji,
    flux_storage_1,
    flux_storage_2,
    flux
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

# @gif for i in eachindex(sol.t)
#     plot(x, getindex.(sol.u[i], 1), leg=false, ylims=(.5, 1.5))
#     plot!(x, getindex.(initial_condition.(x, sol.t[i]), 1), leg=false, ylims=(.5, 1.5))
# end

cons = prim2cons.(sol.u[end], equations)
v1 = getindex.(cons, 2)
v2 = getindex.(cons, 3)

heatmap(reshape(compute_vorticity(v1, v2), num_nodes, num_nodes)',
    c=:inferno,
    aspect_ratio=:equal,
    clim=(-30, 30))