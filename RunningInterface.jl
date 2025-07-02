using Trixi
using SummationByPartsOperators
using OrdinaryDiffEq
using Plots
using LinearAlgebra
using SparseArrays
using FFTW
include("initial_conditions.jl")
include("L2_knapsack.jl")
include("L2_knapsack_maximizer.jl")

accuracy_order = 4
num_nodes = 1024

timestepper = RK4()
abstol = 1e-6
reltol = 1e-4
dt = 1e-4
adaptive = true
saveat = 1e-2

entropy_inequality = :semi_local # which inequality will we enforce
blend = :viscosity # how to blend together schemes to satisfy the inequality

filter_strength = 0.
# 5e-5 for advection, 3, 1000 with buzz
# 5e-4 for burgers, 3, 1000
# 5e-6 for shu osh, 6, 1024

volume_flux = flux_central
low_order_volume_flux = flux_hllc

preserve_positivity = -1

knapsack_shock_capturing = -1

knapsack = QuadraticKnapsackMinimizer{Float64}