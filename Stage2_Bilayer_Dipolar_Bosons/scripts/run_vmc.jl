# Stage2_Bilayer_Dipolar_Bosons/scripts/run_vmc.jl
#
# Expects from main_VMC.jl: num_part, nr0_sq, h, L, r_min, Δ_shoot, tol_shoot,
# R0_opt, num_steps_production, stage_dir (stage1_dir already used by
# optimize_R0.jl, run just before this)

println("Optimal R0 = ", R0_opt)
println("Number of steps: ", num_steps_production)

N_half = num_part ÷ 2
R_match = load_Rmatch(stage1_dir, N_half, L)
Constants = calculate_constants(L, R_match)

# Rebuild the interpolants at R0_opt — optimize_R0 only returned the
# summary numbers (R0_opt, E_opt, err_opt, energy_b_opt), not itp_u/
# itp_up/itp_upp themselves, so this is a genuine second build_fAB call,
# not a redundant recomputation of something already available.
fAB_result = build_fAB(h, R0_opt, r_min, Δ_shoot, tol_shoot, nr0_sq, num_part)
fAB_result === nothing && error("build_fAB failed at R0_opt = $R0_opt — cannot run production VMC.")
_, _, itp_u, itp_up, itp_upp, energy_b = fAB_result

trial = BilayerTrial(; R_match, Constants, R0=R0_opt, itp_u, itp_up, itp_upp, h)

println("\n--- Tuning delta ---")
coords = init_random_config((A=N_half, B=N_half), L)
delta, coords = tune_delta(trial, coords, L; target_ratio=0.5, num_tune_steps=10^4, block_size=100)
println("Tuned delta = ", delta)

println("\n--- Running production VMC ($num_steps_production steps) ---")
result = metropolis(trial, (A=N_half, B=N_half), num_steps_production, delta, L;
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

plateau_std       = detect_plateau(block_sizes, sigmas, window_size=4, rtol=0.05)
plateau_drift     = detect_plateau(block_sizes, sigmas_drift, window_size=4, rtol=0.05)
plateau_laplacian = detect_plateau(block_sizes, sigmas_laplacian, window_size=4, rtol=0.05)

avg_energy, sigma                     = blocking_statistics(result.energies, plateau_std)
avg_energy_drift, sigma_drift         = blocking_statistics(result.energies_drift, plateau_drift)
avg_energy_laplacian, sigma_laplacian = blocking_statistics(result.energies_laplacian, plateau_laplacian)

tail = tail_energy(nr0_sq, num_part, h)
println("\n--- Results ---")
println("E/N = ", avg_energy/num_part, " ± ", sigma/num_part)
println("E/N + tail = ", avg_energy/num_part + tail, " ± ", sigma/num_part)
println("energy_b (dimer binding energy) = ", energy_b)
println("Acceptance ratio = ", result.acceptance_ratio)

path = result_path(stage_dir, "VMC", (N=num_part, nr0sq=nr0_sq, h=h))
save_run(path, (; num_part, nr0_sq, h, L, R_match, R0=R0_opt, energy_b, num_steps_production),
         (; result..., avg_energy, sigma, avg_energy_drift, sigma_drift,
            avg_energy_laplacian, sigma_laplacian, plateau_std, plateau_drift, plateau_laplacian,
            block_sizes, sigmas, sigmas_drift, sigmas_laplacian))
println("Saved results to: ", path)

if plot_vmc_convergence_diagnostic
    p_conv = plot_convergence(block_sizes, sigmas, sigmas_drift, sigmas_laplacian,
                               plateau_std, plateau_drift, plateau_laplacian;
                               title="VMC convergence, N=$(num_part), nr0²=$(nr0_sq), h=$(h)")
    display(p_conv)
    savefig(p_conv, joinpath(stage_dir, "data", "plots", "vmc_convergence_N$(num_part)_nr0sq$(nr0_sq)_h$(h).pdf"))
end

