using Trixi
using SummationByPartsOperators
using OrdinaryDiffEq
using Plots
using LinearAlgebra
using SparseArrays
using FFTW
include("initial_conditions.jl")
include("../L2_knapsack.jl")
include("../L2_knapsack_maximizer.jl")

accuracy_order = 4
num_nodes = 256

timestepper = SSPRK43()
abstol = 1e-6
reltol = 1e-4
dt = 1e-4
adaptive = true
saveat = 1.

entropy_inequality = :semi_local # which inequality will we enforce
blend = :viscosity # how to blend together schemes to satisfy the inequality

low_order_volume_flux = flux_hllc

preserve_positivity = -1

knapsack = QuadraticKnapsackMinimizer{Float64}