# Stage1_2D_Dipol_System/config.jl
#
# Single source of truth for Stage1's run parameters. Included by
# main_VMC.jl, main_DMC.jl, and final_plot.jl, so all three always agree
# on num_part/nr0_sq_vals — change the list once here, not in three places.

num_part    = 30
nr0_sq_vals = [0.5, 16.0, 32.0, 48.0, 64.0, 96.0, 128.0, 196.0, 256.0]
stage_dir   = @__DIR__

plot_Rmatch_sweep_diagnostic = false
plot_dmc_trace_diagnostic    = false