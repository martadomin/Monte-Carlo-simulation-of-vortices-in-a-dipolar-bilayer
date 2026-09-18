# Stage3_Vortex_Excitations/scripts/run_vmc.jl
#
# Production VMC for one vortex offset d — no optimization, R_match/R0
# loaded, vortex core positions built from d (A fixed at box center,
# B offset by d along x). Expects from main_VMC.jl: num_part, nr0_sq, h,
# lA, lB, d, L, r_min, Δ_shoot, tol_shoot, num_steps_production, stage_dir,
# stage1_dir, stage2_dir

println("--- d = $d ---")

N_half = num_part ÷ 2
R_match = load_Rmatch(stage1_dir, N_half, L)
Constants = calculate_constants(L, R_match)
R0 = load_R0(stage2_dir, num_part, nr0_sq, h)

fAB_result = build_fAB(h, R0, r_min, Δ_shoot, tol_shoot, nr0_sq, num_part)
fAB_result === nothing && error("build_fAB failed at R0=$R0 — cannot run production VMC.")
_, _, itp_u, itp_up, itp_upp, _ = fAB_result

x_vortex_A, y_vortex_A = L/2, L/2
x_vortex_B, y_vortex_B = L/2 + d, L/2

trial = VortexBilayerTrial(; R_match, Constants, R0, itp_u, itp_up, itp_upp, h, lA, lB,
                              x_vortex_A, y_vortex_A, x_vortex_B, y_vortex_B)

coords = init_random_config((A=N_half, B=N_half), L)
delta, coords = tune_delta(trial, coords, L; target_ratio=0.5, num_tune_steps=10^4, block_size=100)
println("Tuned delta = ", delta)

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

println("E/N = ", avg_energy/num_part, " ± ", sigma/num_part)
println("Acceptance ratio = ", result.acceptance_ratio)

if plot_vmc_convergence_diagnostic
    p_conv = plot_convergence(block_sizes, sigmas, sigmas_drift, sigmas_laplacian,
                               plateau_std, plateau_drift, plateau_laplacian;
                               title="VMC convergence, N=$(num_part), nr0²=$(nr0_sq), h=$(h), d=$(round(d,digits=4))")
    display(p_conv)
    savefig(p_conv, joinpath(stage_dir, "data", "plots",
            "vmc_convergence_N$(num_part)_nr0sq$(nr0_sq)_h$(h)_d$(round(d,digits=4)).pdf"))
end

path = result_path(stage_dir, "VMC_vortex", (N=num_part, nr0sq=nr0_sq, h=h, lA=lA, lB=lB, d=d))
save_run(path, (; num_part, nr0_sq, h, lA, lB, d, L, R_match, R0, num_steps_production),
         (; result..., avg_energy, sigma, avg_energy_drift, sigma_drift,
            avg_energy_laplacian, sigma_laplacian, plateau_std, plateau_drift, plateau_laplacian,
            block_sizes, sigmas, sigmas_drift, sigmas_laplacian))
println("Saved to: ", path)