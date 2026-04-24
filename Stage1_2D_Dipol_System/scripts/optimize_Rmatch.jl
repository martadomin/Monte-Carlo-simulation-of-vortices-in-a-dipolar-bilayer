using Plots, LaTeXStrings

#First rough sweep to find the approximate location of the optimal R_match:
R_match_vals_coarse = collect(LinRange(0.1*L/2, 0.9*L/2, 20))
println("\n--- Coarse sweep ($(length(R_match_vals_coarse)) points, $num_steps_coarse steps) ---")
results_coarse = sweep_Rmatch(L, num_part, nr0_sq, 
                               R_match_vals_coarse, num_steps_coarse)
R_opt_rough = results_coarse.R_match_vals[argmin(results_coarse.energies)]
println("Rough optimal R_match = $R_opt_rough")

#Second fine sweep around the optimal R_match found in the coarse sweep:
R_match_vals_fine = collect(LinRange(0.7*R_opt_rough, 1.3*R_opt_rough, 20))
println("\n--- Fine sweep ($(length(R_match_vals_fine)) points, $num_steps_fine steps) ---")
results_fine = sweep_Rmatch(L, num_part, nr0_sq, 
                             R_match_vals_fine, num_steps_fine)
R_opt = results_fine.R_match_vals[argmin(results_fine.energies)]

println("\n--- Results ---")
println("Optimal R_match = $R_opt")
println("Optimal energy per particle = ", minimum(results_fine.energies))

# Save results to file
results_path = joinpath(@__DIR__, "..", "data", "sweep_results",
               "Rmatch_sweep_N$(num_part)_nr0sq$(nr0_sq).txt")
open(results_path, "w") do io
    println(io, "# Coarse sweep")
    println(io, "R_match\tEnergy\tDrift_Energy\tLaplacian_Energy")
    for (r, e, e_drift, e_laplacian) in zip(results_coarse.R_match_vals, results_coarse.energies, results_coarse.energies_drift, results_coarse.energies_laplacian)
        println(io, "$(r)\t$(e)\t$(e_drift)\t$(e_laplacian)")
    end
    println(io, "\n# Fine sweep")
    println(io, "R_match\tEnergy\tDrift_Energy\tLaplacian_Energy")
    for (r, e, e_drift, e_laplacian) in zip(results_fine.R_match_vals, results_fine.energies, results_fine.energies_drift, results_fine.energies_laplacian)
        println(io, "$(r)\t$(e)\t$(e_drift)\t$(e_laplacian)")
    end
    println(io, "\n# Optimal R_match")
    println(io, "R_opt\tE_opt")
    println(io, "$(R_opt)\t$(minimum(results_fine.energies))")
end
println("Saved sweep results to: ", results_path)

R_match = R_opt  # make R_match available globally for run_vmc.jl
println("R_match set to R_opt = $R_match")