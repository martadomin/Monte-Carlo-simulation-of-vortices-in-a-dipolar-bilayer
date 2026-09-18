# Stage1_2D_Dipol_System/scripts/tune_delta_dmc.jl
#
# Expects from main_DMC.jl: num_part, nr0_sq, L, R_opt, stage_dir

Constants = calculate_constants(L, R_opt)
trial = SingleLayerTrial(; R_match=R_opt, Constants=Constants)
coords = init_random_config((A=num_part,), L)

println("\n--- Tuning delta for DMC (with VMC move_all criterion) ---")
delta_opt, Δτ_DMC, acceptance = tune_delta_dmc(trial, coords, L)
println("Optimal delta = $delta_opt  (acceptance = $(round(acceptance*100, digits=2))%)")
println("For the DMC run: Δτ = δ²/(2D) = $Δτ_DMC")

path = result_path(stage_dir, "tune_delta_dmc", (N=num_part, nr0sq=nr0_sq))
save_run(path, (; num_part, nr0_sq, L, R_opt), (; delta_opt, Δτ_DMC, acceptance))
println("Saved to: ", path)