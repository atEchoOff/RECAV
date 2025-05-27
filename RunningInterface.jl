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

blend = :semi_local_entropy_knapsack
blending_strat = :fft

filter_strength = 5e-6
# 5e-5 for advection, 3, 1000 with buzz
# 5e-4 for burgers, 3, 1000
# 5e-6 for shu osh, 6, 1024

volume_flux = flux_central

knapsack = QuadraticKnapsackMinimizer{Float64}