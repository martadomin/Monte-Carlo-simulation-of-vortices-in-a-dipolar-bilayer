# Stage1_2D_Dipol_System/scripts/run_vmc.jl
#
# Expects from main_VMC.jl: num_part, nr0_sq, L, R_opt, num_steps_production,
# stage_dir (already defined by optimize_Rmatch.jl, run just before this)

println("Optimal R_match = ", R_opt)
println("Number of steps: ", num_steps_production)

Constants = calculate_constants(L, R_opt)
trial = SingleLayerTrial(; R_match=R_opt, Constants=Constants)

println("\n--- Tuning delta ---")
coords = init_random_config((A=num_part,), L)
delta, coords = tune_delta(trial, coords, L; target_ratio=0.5, num_tune_steps=10^4, block_size=100)
println("Tuned delta = ", delta)

println("\n--- Running production VMC ($num_steps_production steps) ---")
result = metropolis(trial, (A=num_part,), num_steps_production, delta, L;
                     coords_init=coords, final_energy_plot=false,
                     plot_every=10^3, progress=true, num_bins=100)

block_sizes = [10,20,30,40,50,100,150,200,300,400,500,600,700,800,900,
               1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000]
sigmas, sigmas_drift, sigmas_laplacian = Float64[], Float64[], Float64[]
for B in block_sizes
    _, s  = blocking_statistics(result.energies, B)
    _, sd = blocking_statistics(result.energies_drift, B)
    _, sl = blocking_statistics(result.energies_laplacian, B)
    push!(sigmas, s); push!(sigmas_drift, sd); push!(sigmas_laplacian, sl)
end

println("\n--- Automating Plateau Detection ---")
plateau_std       = detect_plateau(block_sizes, sigmas, window_size=4, rtol=0.05)
plateau_drift     = detect_plateau(block_sizes, sigmas_drift, window_size=4, rtol=0.05)
plateau_laplacian = detect_plateau(block_sizes, sigmas_laplacian, window_size=4, rtol=0.05)
println("Detected plateau block size for standard estimator:  ", plateau_std)
println("Detected plateau block size for drift estimator:     ", plateau_drift)
println("Detected plateau block size for laplacian estimator: ", plateau_laplacian)

avg_energy, sigma                     = blocking_statistics(result.energies, plateau_std)
avg_energy_drift, sigma_drift         = blocking_statistics(result.energies_drift, plateau_drift)
avg_energy_laplacian, sigma_laplacian = blocking_statistics(result.energies_laplacian, plateau_laplacian)

nr0_sq_32 = num_part * nr0_sq^(3/2)
println("\n--- Results ---")
println("E/N/(nr0^2)^(3/2) ± σ = ", avg_energy/nr0_sq_32, " ± ", sigma/nr0_sq_32)
println("E_drift/N/(nr0^2)^(3/2) ± σ = ", avg_energy_drift/nr0_sq_32, " ± ", sigma_drift/nr0_sq_32)
println("E_laplacian/N/(nr0^2)^(3/2) ± σ = ", avg_energy_laplacian/nr0_sq_32, " ± ", sigma_laplacian/nr0_sq_32)
println("Acceptance ratio = ", result.acceptance_ratio)

path = result_path(stage_dir, "VMC", (N=num_part, nr0sq=nr0_sq))
save_run(path, (; num_part, nr0_sq, L, R_opt=R_opt, num_steps_production),
         (; result..., avg_energy, sigma, avg_energy_drift, sigma_drift,
            avg_energy_laplacian, sigma_laplacian, plateau_std, plateau_drift, plateau_laplacian))
println("Saved results to: ", path)

if plot_vmc_convergence_diagnostic
    p_conv = plot_convergence(block_sizes, sigmas, sigmas_drift, sigmas_laplacian,
                               plateau_std, plateau_drift, plateau_laplacian;
                               title="VMC convergence, N=$(num_part), nr0²=$(nr0_sq)")
    display(p_conv)
    savefig(p_conv, joinpath(stage_dir, "data", "plots", "vmc_convergence_N$(num_part)_nr0sq$(nr0_sq).pdf"))
end