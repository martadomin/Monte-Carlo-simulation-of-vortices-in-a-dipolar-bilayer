# Stage3_Vortex_Excitations/config.jl

num_part = 60
nr0_sq   = 0.5
h        = 0.7

lA = 1.0
lB = 1.0    # independent, per your earlier confirmation — set differently if needed

r_min       = 1e-6
Δ_shoot     = 1e-4
tol_shoot   = 1e-10

d_vals = collect(0.0:0.05*(sqrt(num_part/nr0_sq)/4):sqrt(num_part/nr0_sq)/4)
# a single d: d_vals = [0.3]

num_steps_production = 10^6
num_steps_dmc         = 10^5

quadratic = true
num_walkers_dtau_study = 200
num_walkers_vals       = [20, 30, 40, 50, 75, 150, 200, 300, 400]

plot_vmc_convergence_diagnostic = false
plot_dmc_trace_diagnostic       = false

stage_dir  = @__DIR__
stage1_dir = joinpath(@__DIR__, "..", "Stage1_2D_Dipol_System")
stage2_dir = joinpath(@__DIR__, "..", "Stage2_Bilayer_Dipolar_Bosons")