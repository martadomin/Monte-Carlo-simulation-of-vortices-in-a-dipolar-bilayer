# Stage3_Vortex_Excitations/scripts/plot_vortex_observables.jl
#
# Loads matching VMC_vortex and DMC_vortex runs, normalizes both, computes
# the extrapolated (linear + ratio) density and g(r), plots all three
# (VMC, DMC mixed, extrapolated) for comparison.

using Plots

vmc_run = load_run(result_path(stage_dir, "VMC_vortex", (N=num_part, nr0sq=nr0_sq, h=h, lA=lA, lB=lB, d=d)); run=1)
dmc_run = load_run(result_path(stage_dir, "DMC_vortex", (N=num_part, nr0sq=nr0_sq, h=h, lA=lA, lB=lB, d=d)); run=1)
obs_vmc, obs_dmc = vmc_run.result.observables, dmc_run.result.observables
N_half = num_part ÷ 2

xy_bins = normalize_density!(obs_vmc.n_xy_A, obs_vmc.n_samples, L)
normalize_density!(obs_vmc.n_xy_B, obs_vmc.n_samples, L)
normalize_density!(obs_dmc.n_xy_A, obs_dmc.n_samples, L)
normalize_density!(obs_dmc.n_xy_B, obs_dmc.n_samples, L)

r_vals = normalize_gr!(obs_vmc.gAA_r, N_half, L, obs_vmc.n_samples)
normalize_gr!(obs_vmc.gBB_r, N_half, L, obs_vmc.n_samples); normalize_gr!(obs_vmc.gAB_r, N_half, L, obs_vmc.n_samples)
normalize_gr!(obs_dmc.gAA_r, N_half, L, obs_dmc.n_samples)
normalize_gr!(obs_dmc.gBB_r, N_half, L, obs_dmc.n_samples); normalize_gr!(obs_dmc.gAB_r, N_half, L, obs_dmc.n_samples)

# Extrapolated (ratio form, for positivity) — density and g(r)
n_A_extrap = extrapolated_estimator_ratio(obs_dmc.n_xy_A, obs_vmc.n_xy_A)
n_B_extrap = extrapolated_estimator_ratio(obs_dmc.n_xy_B, obs_vmc.n_xy_B)
gAA_extrap = extrapolated_estimator_ratio(obs_dmc.gAA_r, obs_vmc.gAA_r)
gBB_extrap = extrapolated_estimator_ratio(obs_dmc.gBB_r, obs_vmc.gBB_r)
gAB_extrap = extrapolated_estimator_ratio(obs_dmc.gAB_r, obs_vmc.gAB_r)

p_nA = plot_density_map(n_A_extrap, xy_bins; title="n_A extrapolated, d=$(round(d,digits=4))")
p_nB = plot_density_map(n_B_extrap, xy_bins; title="n_B extrapolated, d=$(round(d,digits=4))")

p_gr = plot_gr(r_vals, obs_vmc.gAA_r; label="g_AA (VMC)")
plot!(p_gr, r_vals, obs_dmc.gAA_r; label="g_AA (DMC mixed)")
plot!(p_gr, r_vals, gAA_extrap; label="g_AA (extrapolated)", linestyle=:dash)

mkpath(joinpath(stage_dir, "data", "plots"))

display(p_nA); display(p_nB); display(p_gr)
savefig(p_nA, joinpath(stage_dir, "data", "plots", "density_A_extrap_N$(num_part)_nr0sq$(nr0_sq)_h$(h)_d$(round(d,digits=4)).pdf"))
savefig(p_nB, joinpath(stage_dir, "data", "plots", "density_B_extrap_N$(num_part)_nr0sq$(nr0_sq)_h$(h)_d$(round(d,digits=4)).pdf"))
savefig(p_gr, joinpath(stage_dir, "data", "plots", "gr_extrap_N$(num_part)_nr0sq$(nr0_sq)_h$(h)_d$(round(d,digits=4)).pdf"))
println("Saved extrapolated observable plots for d=$d")