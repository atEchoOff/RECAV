using Trixi
using SummationByPartsOperators
using OrdinaryDiffEq
using Plots
using LinearAlgebra
include("initial_conditions.jl")
include("L2_knapsack.jl")

accuracy_order = 3
num_nodes = 1000

timestepper = RK4()
abstol = 1e-6
reltol = 1e-4
dt = 1e-4
adaptive = true
saveat = 1e-2

blend = :loworder
blending_strat = :nodal

filter_strength = 0.
# 5e-5 for advection, 3, 1000 with buzz
# 5e-4 for burgers, 3, 1000

volume_flux = flux_central

knapsack = QuadraticKnapsackMinimizer{Float64}