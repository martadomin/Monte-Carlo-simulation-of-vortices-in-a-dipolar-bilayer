# scripts/run_vmc.jl
# Expects from main.jl: num_part, nr0_sq, L, R_opt, num_steps_production

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
energies,
energies_drift,
E_tot,
E_sq,
E_tot_drift
E_kin,
E_int
acceptance_ratio = metropolis(num_part, num_steps_production, delta, L, R_opt, Constants;
                                            x_init=x_init, y_init=y_init,
                                            final_energy_plot=false,
                                            plot_every=10^3,

                                            progress=true)
                                            
# Block averaging to get final energy estimates
block_sizes = [10, 20, 30, 40, 50, 100, 150, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1500, 1600, 1700, 1800, 1900, 2000]
sigmas = Float64[]
sigmas_drift = Float64[]
sigmas_laplacian = Float64[]
for B in block_sizes
    _, sigma = blocking_statistics(energies_vmc, B)
    _, sigma_drift = blocking_statistics(energies_drift_vmc, B)
    _, sigma_laplacian = blocking_statistics(energies_laplacian_vmc, B)
    push!(sigmas, sigma)
    push!(sigmas_drift, sigma_drift)
    push!(sigmas_laplacian, sigma_laplacian)
end

# Show plot and ask for plateau block size
p = plot(block_sizes, sigmas,
        marker=:circle,
        xlabel="Block size",
        ylabel="Standard error",
        title="Blocking analysis, R_match = $(round(R_match, digits=3))",
        linewidth=2,
        xticks = block_sizes,
        xrotation=45)

plot!(block_sizes, sigmas_drift,
    marker=:square,
    linewidth=2)
display(p)

println("\nR_match = $(round(R_match, digits=3))")
println("Enter plateau block size for standard estimator: ")
plateau_std = parse(Int, readline())
println("Enter plateau block size for drift estimator: ")
plateau_drift = parse(Int, readline())
println("Enter plateau block size for laplacian estimator: ")
plateau_laplacian = parse(Int, readline())
# Get error at chosen block size
avg_energy, sigma = blocking_statistics(energies_vmc, plateau_std)
avg_energy_drift, sigma_drift = blocking_statistics(energies_drift_vmc, plateau_drift)
avg_energy_laplacian, sigma_laplacian = blocking_statistics(energies_laplacian_vmc, plateau_laplacian)

nr0_sq_32 = num_part * nr0_sq^(-3/2)

E_tot = avg_energy/nr0_sq_32
E_drift = avg_energy_drift/nr0_sq_32
E_laplacian = avg_energy_laplacian/nr0_sq_32

println("\n--- Results ---")
println("E/N/(nr0^2)^(3/2)     = ", E_tot/nr0_sq_32)
println("E_drift/N/(nr0^2)^(3/2) = ", E_drift/nr0_sq_32)
println("E_laplacian/N/(nr0^2)^(3/2) = ", E_laplacian/nr0_sq_32)
println("E_kin/N/(nr0^2)^(3/2) = ", E_kin/nr0_sq_32)
println("E_int/N/(nr0^2)^(3/2) = ", E_int/nr0_sq_32)
println("Acceptance ratio       = ", acceptance_ratio)

# Save results
results_path = joinpath(@__DIR__, "..", "data", "results",
               "vmc_N$(num_part)_nr0sq$(nr0_sq).txt")
open(results_path, "w") do io
    println(io, "num_part\tnr0_sq\tL\tR_opt\tE_tot\tE_kin\tE_int\tAcceptance")
    println(io, "$(num_part)\t$(nr0_sq)\t$(L)\t$(R_opt)\t$(E_tot/nr0_sq_32)\t$(E_kin/nr0_sq_32)\t$(E_int/nr0_sq_32)\t$(acceptance_ratio)")
end
println("Saved results to: ", results_path)