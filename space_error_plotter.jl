using LsqFit
using LaTeXStrings

function fit_to_order(x_data, y_data, order)
    @. model(x, p) = p[1] * x ^ order

    return curve_fit(model, x_data, y_data, [.5]).param[1]
end

scatter(Ki2, Li2, xaxis=:log, yaxis=:log, label = L"$N$ = 2, est. order = " * "$(Polynomials.fit(log.(Ki2), log.(Li2), 1)[1])")
scatter!(Ki3, Li3, label = L"$N$ = 3, est. order = " * "$(Polynomials.fit(log.(Ki3), log.(Li3), 1)[1])")
scatter!(Ki4, Li4, label = L"$N$ = 4, est. order = " * "$(Polynomials.fit(log.(Ki4), log.(Li4), 1)[1])")
scatter!(Ki5, Li5, label = L"$N$ = 5, est. order = " * "$(Polynomials.fit(log.(Ki5), log.(Li5), 1)[1])")
scatter!(Ki6, Li6, label = L"$N$ = 6, est. order = " * "$(Polynomials.fit(log.(Ki6), log.(Li6), 1)[1])")

plot!(legend=:topleft)

xlabel!(L"\Delta x")
ylabel!(L"$L^2$ Error")

o2 = fit_to_order(Ki2, Li2, 2)
o3 = fit_to_order(Ki3, Li3, 4)
o4 = fit_to_order(Ki4, Li4, 4)
o5 = fit_to_order(Ki5, Li5, 6)
o6 = fit_to_order(Ki6, Li6, 6)

plot!()

# plot!(annotations=([2e-2], [1e-3], text(L"$\mathcal{O}\left(\Delta x^2\right)$", :red)))
# plot!(Ki2, o2 * (Ki2 .^ 2), label = "", c=:red, lw = 2)

# # plot!(annotations=([5e-2], [2e-5], text(L"$\mathcal{O}\left(\Delta x^3\right)$", :red)))
# plot!(Ki3, o3 * (Ki3 .^ 4), label = "", c=:red, lw = 2)

# # plot!(annotations=([1.4e-1], [5e-6], text(L"$\mathcal{O}\left(\Delta x^4\right)$", :red)))
# plot!(Ki4, o4 * (Ki4 .^ 4), label = "", c=:red, lw = 2)

# # plot!(annotations=([2.3e-1], [1e-6], text(L"$\mathcal{O}\left(\Delta x^5\right)$", :red)))
# plot!(Ki5, o5 * (Ki5 .^ 6), label = "", c=:red, lw = 2)

# plot!(Ki6, o6 * (Ki6 .^ 6), label = "", c=:red, lw = 2)