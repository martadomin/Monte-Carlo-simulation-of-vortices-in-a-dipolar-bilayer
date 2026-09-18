# Stage1_2D_Dipol_System/scripts/final_plot_DMC.jl
#
# Loads every DMC result saved by num_walkers_convergence.jl for this
# num_part, plots E/N vs nr0_sq against the Astrakharchik 2007 fit.
# Expects: num_part, nr0_sq_vals, stage_dir.

using Plots, LaTeXStrings

E_dmc, error_E_dmc = Float64[], Float64[]
for nr0_sq in nr0_sq_vals
    path = result_path(stage_dir, "DMC", (N=num_part, nr0sq=nr0_sq))
    run = load_run(path; run=1)
    push!(E_dmc, run.result.E_extrapolated / (num_part * nr0_sq^(3/2)))
    push!(error_E_dmc, run.result.E_extrapolated_err / (num_part * nr0_sq^(3/2)))
end

a1, a2, a3 = 4.536, 4.38, 1.2
E_tail = tail_energy(num_part)
nr0_range = collect(LinRange(10.0, 300.0, 500))
E_paper_normalized = (@. a1*nr0_range^(3/2) + a2*nr0_range^(5/4) + a3*nr0_range^(1/2)) ./ nr0_range.^(3/2)

gr()
p = plot(nr0_range, E_paper_normalized,
          label=L"\mathrm{DMC\ fit\ (Astrakharchik\ 2007)}",
          xlabel=L"nr_0^2", ylabel=L"E/N \cdot (nr_0^2)^{-3/2}",
          title=L"\mathrm{DMC},\ N = %$(num_part)",
          linewidth=1, color=:black, framestyle=:box, legend=:topright)
scatter!(nr0_sq_vals, E_dmc .+ E_tail, yerror=error_E_dmc,
         label=L"\mathrm{DMC\ results,\ }N = %$(num_part)", marker=:square, markersize=6, color=:green)

savefig(joinpath(stage_dir, "data", "plots", "plot_dmc_N$(num_part).pdf"))
println("Saved DMC plot")
display(p)