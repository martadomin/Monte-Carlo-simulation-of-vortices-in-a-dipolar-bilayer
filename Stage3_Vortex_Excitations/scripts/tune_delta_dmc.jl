# Stage3_Vortex_Excitations/scripts/tune_delta_dmc.jl
#
# Expects from main_DMC.jl: num_part, nr0_sq, h, lA, lB, d, L, r_min,
# Δ_shoot, tol_shoot, stage_dir, stage1_dir, stage2_dir

N_half = num_part ÷ 2
R_match = load_Rmatch(stage1_dir, N_half, L)
Constants = calculate_constants(L, R_match)
R0 = load_R0(stage2_dir, num_part, nr0_sq, h)

fAB_result = build_fAB(h, R0, r_min, Δ_shoot, tol_shoot, nr0_sq, num_part)
fAB_result === nothing && error("build_fAB failed at R0=$R0 — cannot tune DMC delta.")
_, _, itp_u, itp_up, itp_upp, _ = fAB_result

x_vortex_A, y_vortex_A = L/2, L/2
x_vortex_B, y_vortex_B = L/2 + d, L/2

trial = VortexBilayerTrial(; R_match, Constants, R0, itp_u, itp_up, itp_upp, h, lA, lB,
                              x_vortex_A, y_vortex_A, x_vortex_B, y_vortex_B)
coords = init_random_config((A=N_half, B=N_half), L)

println("\n--- Tuning delta for DMC (move_all), d=$d ---")
delta_opt, Δτ_DMC, acceptance = tune_delta_dmc(trial, coords, L)
println("Optimal delta = $delta_opt  (acceptance = $(round(acceptance*100, digits=2))%)")
println("For the DMC run: Δτ = δ²/(2D) = $Δτ_DMC")

path = result_path(stage_dir, "tune_delta_dmc", (N=num_part, nr0sq=nr0_sq, h=h, lA=lA, lB=lB, d=d))
save_run(path, (; num_part, nr0_sq, h, lA, lB, d, L, R_match, R0), (; delta_opt, Δτ_DMC, acceptance))
println("Saved to: ", path)