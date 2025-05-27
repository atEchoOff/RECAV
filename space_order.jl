using Polynomials

MODULE = "density_wave.jl"

accuracy_order = 5
adaptive = false
dt = 1e-4

Li = Float64[]
Ki = Float64[]

i = 16
while length(Li) == 0 || Li[end] > 1e-6
    global num_nodes, i
    num_nodes = floor(Int, i)
    println("Testing K = $num_nodes")
    i *= 2^.2

    push!(Ki, 2 / num_nodes)

    include(MODULE)

    nm = L2_error
    println("Obtained norm $nm")
    push!(Li, nm)
end


println(Polynomials.fit(log.(Ki), log.(Li), 1))
display(Li)