# Stage1_2D_Dipol_System/scripts/optimize_Rmatch.jl
#
# Expects from main_VMC.jl: L, num_part, num_steps_coarse, num_steps_fine

stage_dir = joinpath(@__DIR__, "..")
R_opt = optimize_Rmatch(L, num_part, num_steps_coarse, num_steps_fine, stage_dir)
println("R_match set to R_opt = $R_opt")
