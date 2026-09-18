# Stage2_Bilayer_Dipolar_Bosons/scripts/tune_delta_dmc.jl
#
# Expects from main_DMC.jl: num_part, nr0_sq, h, L, r_min, Δ_shoot, tol_shoot,
# stage_dir, stage1_dir

N_half = num_part ÷ 2
R_match = load_Rmatch(stage1_dir, N_half, L)
Constants = calculate_constants(L, R_match)

# R0 was already optimized by main_VMC.jl for this (nr0_sq, h) — load it
# rather than re-running optimize_R0.
R0_run = load_run(result_path(stage_dir, "R0_optimum", (N=num_part, nr0sq=nr0_sq, h=h)); run=1)
R0_opt = R0_run.result.R0_opt
println("Loaded R0_opt = $R0_opt from Stage2's VMC optimization")

fAB_result = build_fAB(h, R0_opt, r_min, Δ_shoot, tol_shoot, nr0_sq, num_part)
fAB_result === nothing && error("build_fAB failed at R0_opt = $R0_opt — cannot tune DMC delta.")
_, _, itp_u, itp_up, itp_upp, _ = fAB_result

trial = BilayerTrial(; R_match, Constants, R0=R0_opt, itp_u, itp_up, itp_upp, h)
coords = init_random_config((A=N_half, B=N_half), L)

println("\n--- Tuning delta for DMC (move_all) ---")
delta_opt, Δτ_DMC, acceptance = tune_delta_dmc(trial, coords, L)
println("Optimal delta = $delta_opt  (acceptance = $(round(acceptance*100, digits=2))%)")
println("For the DMC run: Δτ = δ²/(2D) = $Δτ_DMC")

path = result_path(stage_dir, "tune_delta_dmc", (N=num_part, nr0sq=nr0_sq, h=h))
save_run(path, (; num_part, nr0_sq, h, L, R_match, R0=R0_opt), (; delta_opt, Δτ_DMC, acceptance))
println("Saved to: ", path)