# Stage1_2D_Dipol_System/scripts/final_plot.jl
#
# Standalone — run directly, does not need main_VMC.jl/main_DMC.jl to have
# run in the same session. Requires both VMC and DMC pipelines to have
# already been run (and saved) for every nr0_sq listed below.

include(joinpath(@__DIR__, "..", "common", "src", "CommonCore.jl"))
include(joinpath(@__DIR__, "src", "wavefunction.jl"))
include(joinpath(@__DIR__, "config.jl"))
using Plots, LaTeXStrings

E_vmc, error_E_vmc = Float64[], Float64[]
E_dmc, error_E_dmc = Float64[], Float64[]
for nr0_sq in nr0_sq_vals
    vmc_run = load_run(result_path(stage_dir, "VMC", (N=num_part, nr0sq=nr0_sq)); run=1)
    push!(E_vmc, vmc_run.result.avg_energy / (num_part * nr0_sq^(3/2)))
    push!(error_E_vmc, vmc_run.result.sigma / (num_part * nr0_sq^(3/2)))

    dmc_run = load_run(result_path(stage_dir, "DMC", (N=num_part, nr0sq=nr0_sq)); run=1)
    push!(E_dmc, dmc_run.result.E_extrapolated / (num_part * nr0_sq^(3/2)))
    push!(error_E_dmc, dmc_run.result.E_extrapolated_err / (num_part * nr0_sq^(3/2)))
end

a1, a2, a3 = 4.536, 4.38, 1.2
E_tail = tail_energy(num_part)
nr0_range = collect(LinRange(10.0, 300.0, 500))
E_paper_normalized = (@. a1*nr0_range^(3/2) + a2*nr0_range^(5/4) + a3*nr0_range^(1/2)) ./ nr0_range.^(3/2)

gr()
p = plot(nr0_range, E_paper_normalized,
          label=L"\mathrm{DMC\ fit\ (Astrakharchik\ 2007)}",
          xlabel=L"nr_0^2", ylabel=L"E/N \cdot (nr_0^2)^{-3/2}",
          title=L"\mathrm{VMC\ vs\ DMC},\ N = %$(num_part)",
          linewidth=1, color=:black, framestyle=:box, legend=:topright)
scatter!(nr0_sq_vals, E_vmc .+ E_tail, yerror=error_E_vmc,
         label=L"\mathrm{VMC,\ }N = %$(num_part)", marker=:circle, markersize=6, color=:red)
scatter!(nr0_sq_vals, E_dmc .+ E_tail, yerror=error_E_dmc,
         label=L"\mathrm{DMC\ results,\ }N = %$(num_part)", marker=:square, markersize=6, color=:green)

savefig(joinpath(stage_dir, "data", "plots", "plot_vmc_vs_dmc_N$(num_part).pdf"))
println("Saved combined plot")
display(p)