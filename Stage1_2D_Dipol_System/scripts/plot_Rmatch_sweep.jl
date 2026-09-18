# Stage1_2D_Dipol_System/scripts/plot_Rmatch_sweep.jl
#
# Plots the fine R_match sweep — E/N vs R_match, all three estimators
# (std/drift/laplacian) overlaid, R_opt marked. Diagnostic of the
# optimization process. Expects: num_part, L, stage_dir

using Plots, LaTeXStrings

run = load_run(result_path(stage_dir, "Rmatch_optimum", (N=num_part, L=L)); run=1)
r = run.result.results_fine   # coarse sweep also available via run.result.results_coarse

p = plot(r.R_match_vals, r.energies ./ num_part; yerror=r.error ./ num_part, label="std",
         marker=:circle, color=:blue, xlabel=L"R_{\mathrm{match}}", ylabel=L"E/N",
         title="R_match sweep, N=$(num_part), L=$(round(L,digits=3))", framestyle=:box)
plot!(r.R_match_vals, r.energies_drift ./ num_part; yerror=r.error_drift ./ num_part,
      label="drift", marker=:square, color=:green)
plot!(r.R_match_vals, r.energies_laplacian ./ num_part; yerror=r.error_laplacian ./ num_part,
      label="laplacian", marker=:diamond, color=:orange)
vline!([run.result.R_opt]; label="R_opt", linestyle=:dash, color=:red)

display(p)
savefig(p, joinpath(stage_dir, "data", "plots", "Rmatch_sweep_N$(num_part)_L$(round(L,digits=3)).pdf"))
println("Saved R_match sweep plot")