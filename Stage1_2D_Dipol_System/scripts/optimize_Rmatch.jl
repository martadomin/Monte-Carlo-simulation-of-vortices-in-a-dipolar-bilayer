include(normpath(joinpath(@__DIR__, "..", "src", "optimization.jl")))
using Plots, LaTeXStrings

#First rough sweep to find the approximate location of the optimal R_match:
R_match_vals_coarse = collect(LinRange(0.1*L/2, 0.9*L/2, 20))
results_coarse = sweep_Rmatch(L, num_part, nr0_sq, 
                               R_match_vals_coarse, num_steps_coarse)
R_opt_rough = results_coarse.R_match_vals[argmin(results_coarse.energies)]

#Second fine sweep around the optimal R_match found in the coarse sweep:
R_match_vals_fine = collect(LinRange(0.7*R_opt_rough, 1.3*R_opt_rough, 50))
results_fine = sweep_Rmatch(L, num_part, nr0_sq, 
                             R_match_vals_fine, num_steps_fine)
R_opt = results_fine.R_match_vals[argmin(results_fine.energies)]

println("Optimal R_match = $R_opt")
println("Optimal energy per particle = ", minimum(results_fine.energies))

# Save results to file
results_path = joinpath(@__DIR__, "..", "data", "sweep_results",
               "Rmatch_sweep_N$(num_part)_nr0sq$(nr0_sq).txt")
open(results_path, "w") do io
    println(io, "# Coarse sweep")
    println(io, "R_match\tEnergy")
    for (r, e) in zip(results_coarse.R_match_vals, results_coarse.energies)
        println(io, "$(r)\t$(e)")
    end
    println(io, "\n# Fine sweep")
    println(io, "R_match\tEnergy")
    for (r, e) in zip(results_fine.R_match_vals, results_fine.energies)
        println(io, "$(r)\t$(e)")
    end
    println(io, "\n# Optimal R_match")
    println(io, "R_opt\tE_opt")
    println(io, "$(R_opt)\t$(minimum(results_fine.energies))")
end
println("Saved sweep results to: ", results_path)

# Plot
gr()
plot(results_coarse.R_match_vals, results_coarse.energies,
     label=L"Coarse sweep",
     xlabel=L"R_{\mathrm{match}}",
     ylabel=L"E/N \cdot (nr_0^2)^{-3/2}",
     title="VMC energy vs \$R_{\\mathrm{match}}\$, \$nr_0^2\$ = $(nr0_sq), N = $(num_part)",
     marker=:circle,
     linewidth=2)
plot!(results_fine.R_match_vals, results_fine.energies,
      label=L"Fine sweep",
      marker=:circle,
      linewidth=2)
vline!([R_opt], linestyle=:dash, color=:red,
       label="R_opt = $(round(R_opt, digits=4))")
savefig(joinpath(@__DIR__, "..", "data", "sweep_results",
        "Rmatch_sweep_N$(num_part)_nr0sq$(nr0_sq).pdf"))