total_entropies = Float64[]

for i in eachindex(sol.t)
    total_entropy = sum(M * entropy.(sol.u[i], equations))
    push!(total_entropies, total_entropy)
end