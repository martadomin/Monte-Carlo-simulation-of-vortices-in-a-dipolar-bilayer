# Stage2_Bilayer_Dipolar_Bosons/scripts/optimize_R0.jl
#
# Expects from main_VMC.jl: num_part, nr0_sq, h, L, r_min, Δ_shoot, tol_shoot,
# n_points_sweep, num_steps_coarse, num_steps_fine, stage_dir, stage1_dir

N_half = num_part ÷ 2
R_match = load_Rmatch(stage1_dir, N_half, L)
Constants = calculate_constants(L, R_match)
println("Loaded R_match = $R_match from Stage1 (N_half=$N_half, L=$L)")

