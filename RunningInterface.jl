using Trixi
using SummationByPartsOperators
using OrdinaryDiffEq
using Plots
using LinearAlgebra
using SparseArrays
using FFTW
include("initial_conditions.jl")
include("L2_knapsack.jl")

accuracy_order = 6
num_nodes = 1024

timestepper = RK4()
abstol = 1e-6
reltol = 1e-4
dt = 1e-4
adaptive = true
saveat = 1e-2

entropy_inequality = :semi_local # which inequality will we enforce
blend = :knapsack # how to blend together schemes to satisfy the inequality
entropy_blend = :grouped # how to treat the low/high order entropy fluxes
blending_strat = :nodal # for global knapsack, blending strat. Otherwise, postprocessing step.

filter_strength = 0.
# 5e-5 for advection, 3, 1000 with buzz
# 5e-4 for burgers, 3, 1000
# 5e-6 for shu osh, 6, 1024

volume_flux = flux_central

knapsack = QuadraticKnapsackMinimizer{Float64}