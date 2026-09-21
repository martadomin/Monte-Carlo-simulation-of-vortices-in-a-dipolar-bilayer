# Stage1_2D_Dipol_System/scripts/plot_Rmatch_sweep.jl
#
# Plots BOTH the coarse and fine R_match sweeps — E/N vs R_match, all
# three estimators overlaid per sweep, R_opt marked. Two panels side by
# side so you can see whether the fine sweep actually zoomed in on the
# coarse sweep's true minimum, or drifted elsewhere due to sweep noise.
# Expects: num_part, L, stage_dir

using Plots, LaTeXStrings

run = load_run(result_path(stage_dir, "Rmatch_optimum", (N=num_part, L=L)); run=1)
coarse, fine = run.result.results_coarse, run.result.results_fine

# plot_Rmatch_sweep.jl — add right after loading, before plotting
println("Coarse energies: ", coarse.energies)
println("Fine energies:   ", fine.energies)

function _sweep_panel(r, title_str)
    p = plot(r.R_match_vals, r.energies ./ num_part; yerror=r.error ./ num_part,
             label="std", marker=:circle, color=:blue, xlabel=L"R_{\mathrm{match}}",
             ylabel=L"E/N", title=title_str, framestyle=:box)
    plot!(r.R_match_vals, r.energies_drift ./ num_part; yerror=r.error_drift ./ num_part,
          label="drift", marker=:square, color=:green)
    plot!(r.R_match_vals, r.energies_laplacian ./ num_part; yerror=r.error_laplacian ./ num_part,
          label="laplacian", marker=:diamond, color=:orange)
    return p
end

p_coarse = _sweep_panel(coarse, "Coarse sweep")
vline!(p_coarse, [run.result.R_opt_rough]; label="R_opt_rough", linestyle=:dash, color=:red)

p_fine = _sweep_panel(fine, "Fine sweep")
vline!(p_fine, [run.result.R_opt]; label="R_opt", linestyle=:dash, color=:red)

p = plot(p_coarse, p_fine; layout=(1,2), size=(1200,500),
         plot_title="R_match sweep, N=$(num_part), L=$(round(L,digits=3))")

display(p)
savefig(p, joinpath(stage_dir, "data", "plots", "Rmatch_sweep_both_N$(num_part)_L$(round(L,digits=3)).pdf"))
println("Saved R_match sweep plot (coarse + fine)")