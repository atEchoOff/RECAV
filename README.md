# (R)ECAV

This is the testing repo for the paper *Entropy stable finite difference (ESFD) methods via entropy correction artificial viscosity (ECAV) and knapsack limiting (KL) techniques*.

In order to run tests, first run either 1D/RunningInterface.jl or 2D/RunningInterface.jl. Then, run the corresponding problem. Here is an explanation of each tunable setting:

- `accuracy_order` defines the degree of the chosen finite difference stencil. In the case of a smooth problem, it will also correspond to a lower bound of degree of spatial convergence.
- `num_nodes` defines the number of nodes along each axis. In the case of two dimensions, the total number of nodes is the square of `num_nodes`.
- `timestepper` is the chosen timestepper available from `OrdinaryDiffEq.jl`. 
- `abstol` and `reltol` correspond to the absolute and relative tolerances, if `adaptive = true`.
- `dt` corresponds to either the initial timestep if `adaptive = true`, or the fixed timestep if `adaptive = false`.
- `adaptive` is a boolean which enables or disables adaptive timestepping.
- `saveat` defines how often to save a solution snapshot.
- `weak_bcs` toggles between strongly and weakly enforced boundary conditions (only implemented in 1D)
- `entropy_inequality` should either be set to `:semi_local` for the entropy stable schemes, or `:none` otherwise.
- `blend` should be set to `:viscosity` for the ECAV scheme, or `:knapsack` for the KL scheme.
- `potential_blend` should be set to `:relaxed` to use the relaxed entropy inequality (with tau coefficients). It can be set to anything else to use the original entropy inequality. (only implemented in 1D)
- `low_order_volume_flux` sets the flux used for weak boundary conditions (if `weak_bcs = true`) and the low order flux for knapsack limiting (if `blend = :knapsack`).
- `preserve_positivity` enables positivity preservation, if `blend = :knapsack`. If it is set to `-1`, positivity preservation is disabled. If it is between `0` and `1`, it corresponds to the chosen alpha for the relative positivity inequality.
- `knapsack` determines the chosen knapsack limiter. Currently, only `QuadraticKnapsackMinimizer` is implemented.