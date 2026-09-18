# Stage2_Bilayer_Dipolar_Bosons/scripts/plot_R0_sweep.jl
#
# Plots the fine R0 sweep — E/N vs R0, all three estimators (std/drift/
# laplacian) overlaid, R0_opt marked. Diagnostic of the optimization
# process, analogous to Stage1's plot_convergence (std/drift/laplacian
# vs block size). Expects: num_part, nr0_sq, h, stage_dir

using Plots, LaTeXStrings

run = load_run(result_path(stage_dir, "R0_optimum", (N=num_part, nr0sq=nr0_sq, h=h)); run=1)
r = run.result.fine   # coarse sweep also available via run.result.coarse

p = plot(r.R0_vals, r.energies; yerror=r.errors, label="std", marker=:circle, color=:blue,
         xlabel=L"R_0", ylabel=L"E/N", title="R0 sweep, N=$(num_part), nr0²=$(nr0_sq), h=$(h)",
         framestyle=:box)
plot!(r.R0_vals, r.energies_drift; yerror=r.errors_drift, label="drift", marker=:square, color=:green)
plot!(r.R0_vals, r.energies_lap; yerror=r.errors_lap, label="laplacian", marker=:diamond, color=:orange)
vline!([run.result.R0_opt]; label="R0_opt", linestyle=:dash, color=:red)

display(p)
savefig(p, joinpath(stage_dir, "data", "plots", "R0_sweep_N$(num_part)_nr0sq$(nr0_sq)_h$(h).pdf"))
println("Saved R0 sweep plot")