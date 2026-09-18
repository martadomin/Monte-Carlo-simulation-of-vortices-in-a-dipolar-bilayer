# Stage2_Bilayer_Dipolar_Bosons/config.jl
#
# Single source of truth for Stage2's run parameters. Included by
# main_VMC.jl, main_DMC.jl, and any standalone analysis script.

num_part    = 60
nr0_sq_vals = [0.5]
h_vals      = [0.3, 0.4, 0.5, 0.7, 1.0, 1.3, 1.5]

r_min     = 1e-6
Δ_shoot   = 1e-4
tol_shoot = 1e-10

n_points_sweep = 16    # R0 sweep resolution (coarse and fine), matching
                        #   Stage1's own hardcoded 16-point R_match sweeps

stage_dir  = @__DIR__
stage1_dir = joinpath(@__DIR__, "..", "Stage1_2D_Dipol_System")

plot_R0_sweep_diagnostic  = false
plot_dmc_trace_diagnostic = false