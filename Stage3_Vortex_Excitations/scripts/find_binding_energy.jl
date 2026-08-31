# include(normpath(joinpath(@__DIR__, "..", "src", "shooting_method.jl")))
using Plots, LaTeXStrings, Base.Threads

# ================================================================
# MAIN
# ================================================================

R_inf  = L / 2        # large enough that ψ_b has decayed
energy_b_vals = Float64[]

open(normpath(joinpath(@__DIR__, "..", "data", "binding_energy_dimer",
              "dimer_binding_energy.txt")), "w") do io
    println(io, "# h/r0    ε_b")
    for h in h_vals
        Δ_h = min(Δ_shoot, h^(3/2) / 20.0)
        eb  = find_energy_b(r_min, h, R_inf, Δ_h, tol_shoot)
        push!(energy_b_vals, eb)
        # println("h = $h  →  ε_b = $(round(eb, digits=6))")
        println(io, "$(round(h, digits=6))\t$(round(eb, digits=6))")
    end
end

exact_eps_b = Dict(h => eb for (h, eb) in zip(h_vals, energy_b_vals))
# # Plot ε_b(h)
# pgfplotsx()

# p_ε_b = plot(
#     h_vals,
#     abs.(energy_b_vals) ./ 2.0;
#     xlabel = L"h/r_0",
#     ylabel = L"|\epsilon_b|/2",
#     title = "Dimer Binding Energy",
#     legend = false,
#     grid = true,
#     gridalpha = 0.25,
#     framestyle = :box,
#     background_color = :white,
#     foreground_color = :black,
#     size = (900, 600),
#     xlims = (0.18, 0.6),
#     ylims = (0.0, 40),
#     tickfontsize = 11,
#     guidefontsize = 13,
#     titlefontsize = 15,
#     linewidth = 2,
#     color = :dodgerblue,
# )
# display(p_ε_b)
# savefig(p_ε_b, normpath(joinpath(@__DIR__, "..", "data", "plots", "dimer_binding_energy_for_small_h.png")))
# readline()  # Wait for user input before closing the plot
