# scripts/run_vmc.jl
# Expects from main.jl: num_part, nr0_sq, L, R_opt, num_steps_production

using Plots, LaTeXStrings

# Read R_opt from sweep file
rmatch_path = joinpath(@__DIR__, "..", "data", "sweep_results",
              "Rmatch_sweep_N$(num_part)_nr0sq$(nr0_sq).txt")

R_opt_ref = Ref(0.0)
open(rmatch_path, "r") do io
    section = ""
    for line in eachline(io)
        startswith(line, "#") && (section = line; continue)
        isempty(strip(line)) && continue
        line == "R_opt\tE_opt" && continue
        if occursin("Optimal R_match", section)
            vals = parse.(Float64, split(line, "\t"))
            R_opt_ref[] = vals[1]
        end
    end
end
R_opt = R_opt_ref[]
println("Loaded R_opt = $(R_opt) from sweep file")

println("Number of particles: ", num_part)
println("Density nr0^2 = ", nr0_sq)
println("L = ", L)
println("Optimal R_match = ", R_opt)
println("Number of steps: ", num_steps_production)

println("\n--- Tuning delta ---")
Constants = calculate_constants(L, R_opt)
x_coord, y_coord = random_initial_config(num_part, L, "Uniform")
delta, x_init, y_init = tune_delta(x_coord, y_coord, L, R_opt, Constants)
println("Tuned delta = ", delta)

println("\n--- Running production VMC ($num_steps_production steps) ---")
energies_vmc, energies_drift_vmc, energies_laplacian_vmc, E_tot, _, E_drift, E_laplacian, _, _, acceptance_ratio, r_vals, g_r_normalized = metropolis(num_part, 
                                                                                                                            num_steps_production,
                                                                                                                            delta, L, R_opt,
                                                                                                                            Constants;
                                                                                                                            x_init=x_init, y_init=y_init,
                                                                                                                            final_energy_plot=false,
                                                                                                                            plot_every=10^3,
                                                                                                                            progress=true, num_bins=100)
                                            
# Block averaging to get final energy estimates
block_sizes = [10, 20, 30, 40, 50, 100, 150, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000]
sigmas = Float64[]
sigmas_drift = Float64[]
sigmas_laplacian = Float64[]
# for B in block_sizes
#     _, sigma = blocking_statistics(energies_vmc, B)
#     _, sigma_drift = blocking_statistics(energies_drift_vmc, B)
#     _, sigma_laplacian = blocking_statistics(energies_laplacian_vmc, B)
#     push!(sigmas, sigma)
#     push!(sigmas_drift, sigma_drift)
#     push!(sigmas_laplacian, sigma_laplacian)
# end

# p = plot(block_sizes, sigmas,
#         marker=:circle,
#         xlabel="Block size",
#         ylabel="Standard error",
#         title="Blocking analysis, R_opt = $(round(R_opt, digits=4))",
#         linewidth=2,
#         xticks=block_sizes,
#         xrotation=45,
#         label="Standard")
# plot!(block_sizes, sigmas_drift,
#       marker=:square, linewidth=2, label="Drift")
# plot!(block_sizes, sigmas_laplacian,
#       marker=:diamond, linewidth=2, label="Laplacian")

# display(p)

# println("\nR_opt = $(round(R_opt, digits=4))")

# println("Enter plateau block size for standard estimator: ")
# plateau_std = parse(Int, readline())
# println("Enter plateau block size for drift estimator: ")
# plateau_drift = parse(Int, readline())
# println("Enter plateau block size for laplacian estimator: ")
# plateau_laplacian = parse(Int, readline())

plateau_std = 500
plateau_drift = 800
plateau_laplacian = 1000

# Get error at chosen block size
avg_energy, sigma = blocking_statistics(energies_vmc, plateau_std)
avg_energy_drift, sigma_drift = blocking_statistics(energies_drift_vmc, plateau_drift)
avg_energy_laplacian, sigma_laplacian = blocking_statistics(energies_laplacian_vmc, plateau_laplacian)

nr0_sq_32 = num_part * nr0_sq^(3/2)

println("\n--- Results ---")
println("E/N/(nr0^2)^(3/2) ± σ = ", avg_energy/nr0_sq_32, " ± ", sigma/nr0_sq_32)
println("E_drift/N/(nr0^2)^(3/2) ± σ = ", avg_energy_drift/nr0_sq_32 , " ± ", sigma_drift/nr0_sq_32)
println("E_laplacian/N/(nr0^2)^(3/2) ± σ = ", avg_energy_laplacian/nr0_sq_32, " ± ", sigma_laplacian/nr0_sq_32)
println("Acceptance ratio = ", acceptance_ratio)

# Save results
results_path = joinpath(@__DIR__, "..", "data", "results",
               "vmc_N$(num_part)_nr0sq$(nr0_sq).txt")
open(results_path, "w") do io
    println(io, "num_part\tnr0_sq\tL\tR_opt\tE_tot\tError")
    println(io, "$(num_part)\t$(nr0_sq)\t$(L)\t$(Float64(R_opt))\t$(E_tot)\t$(sigma)")
end
println("Saved results to: ", results_path)

# Save g(r)
gr_path = joinpath(@__DIR__, "..", "data", "results",
          "gr_N$(num_part)_nr0sq$(nr0_sq).txt")
open(gr_path, "w") do io
    println(io, "r\tg(r)")
    for (r, g) in zip(r_vals, g_r_normalized)
        println(io, "$(r)\t$(g)")
    end
end
