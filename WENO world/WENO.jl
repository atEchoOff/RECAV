# 1. Import necessary packages
using MethodOfLines
using ModelingToolkit
using DomainSets
using OrdinaryDiffEq

# 2. Define the parameters and variables
# Independent variables
@parameters t x
# Dependent variable
@variables u(..)

# 3. Define the differential operators
Dt = Differential(t)
Dx = Differential(x)

# 4. Specify the 1D linear advection equation
# The advection velocity is set to 1.0
eq = Dt(u(t, x)) + 1.0 * Dx(u(t, x)) ~ 0

# 5. Define the initial and boundary conditions
# The spatial domain is x ∈ [0.0, 1.0].
# For periodic boundary conditions, we state that u at the start of the domain
# equals u at the end of the domain for all time.
# The initial condition is a Gaussian pulse.
bcs = [u(0, x) ~ exp(-10 * x^2), # Made the pulse a bit sharper to see WENO's effect
       u(t, -1.0) ~ u(t, 1.0)]

# 6. Define the spatial and temporal domains
# Spatial domain from 0.0 to 1.0
# Temporal domain from 0.0 to 2.0
domains = [t ∈ Interval(0.0, 1.0),
           x ∈ Interval(-1.0, 1.0)]

# 7. Create the PDE system
@named pdesys = PDESystem(eq, bcs, domains, [t, x], [u(t, x)])

# 8. Define the discretization strategy
# CORRECTED LINE: We now specify the advection scheme to be WENO().
# This uses the 5th-order WENO-JS scheme for the advection term.
discretization = MOLFiniteDifference([x => 1024], t, advection_scheme = WENOScheme())

# 9. Convert the PDE system into an ODE problem
prob = discretize(pdesys, discretization)

# 10. Solve the ODE problem
# We use the Tsit5() solver, which is a good general-purpose solver
sol = solve(prob, RK4(), dt = 1e-4, 
            abstol=1e-6, 
            reltol=1e-4, 
            saveat=1e-2, 
            callback=AliveCallback(alive_interval=1000), 
            adaptive=false)

# 11. Print a success message
println("Solution using 5th-Order WENO successfully computed.")

# To visualize the solution, you can use the Plots.jl package.
# The WENO scheme should preserve the shape of the pulse with very little
# numerical diffusion or oscillation.
#
# First, ensure Plots is installed:
# using Pkg; Pkg.add("Plots")
#
# Then, you can plot the results:
using Plots
#
# # Create a discrete representation of the solution
discrete_x = sol[x]
discrete_t = sol[t]
discrete_u = sol[u(t, x)]
# plot(discrete_x, discrete_u[end, :],
#          ylim=(-0.2, 1.2),
#          xlabel="x", ylabel="u(t,x)",
#          title="Linear Advection with 5th-Order WENO")
#
# # Create an animation to see the periodic behavior
anim = @animate for i in 1:length(discrete_t)
    plot(discrete_x, discrete_u[i, :],
         ylim=(-0.2, 1.2), label="t=$(round(discrete_t[i], digits=2))",
         xlabel="x", ylabel="u(t,x)",
         title="Linear Advection with 5th-Order WENO")
end
#
# # Save the animation as a GIF
gif(anim, "linear_advection_weno.gif", fps = 20)