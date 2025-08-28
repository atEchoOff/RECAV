xmin = -1
xmax = 1
is_periodic = true
T = 4

include("common.jl")

equations = CompressibleEulerEquations2D(1.4)
initial_condition = initial_condition_density_wave_fast

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
    low_order_volume_flux,
    equations, 
    r_H, 
    a, 
    θ, 
    v,
    knapsack_solvers,
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
            callback=AliveCallback(alive_interval=100), 
            adaptive=adaptive)

# @gif for i in eachindex(sol.t)
#     plot(x, getindex.(sol.u[i], 1), leg=false, ylims=(.5, 1.5))
#     plot!(x, getindex.(initial_condition.(x, sol.t[i]), 1), leg=false, ylims=(.5, 1.5))
# end

L2_error = sqrt(sum(M * map(x -> sum(x.^2), initial_condition.(x, y, sol.t[end]) - sol.u[end])))
println(L2_error)

heatmap(reshape(getindex.(sol.u[end], 1), num_nodes, num_nodes)',
    c=:inferno,
    aspect_ratio=:equal)