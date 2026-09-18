# Stage2_Bilayer_Dipolar_Bosons/final_fig1.jl
# Standalone — same requirements as final_plot_inset_combined.jl.
# Fig1 only: E/N - ε_b/2 vs h/r₀, VMC and DMC together.

include(joinpath(@__DIR__, "..", "common", "src", "CommonCore.jl"))
include(joinpath(@__DIR__, "src", "wavefunction.jl"))
include(joinpath(@__DIR__, "src", "h_sweep_analysis.jl"))
include(joinpath(@__DIR__, "config.jl"))
using Plots, LaTeXStrings

nr0_sq = nr0_sq_vals[1]

h_vmc, E_vmc, err_vmc = load_h_sweep(num_part, nr0_sq, h_vals, stage_dir)
h_dmc, E_dmc, err_dmc = load_h_sweep_dmc(num_part, nr0_sq, h_vals, stage_dir)
eb = exact_binding_energy_interp(stage_dir)
E_single = single_layer_reference(num_part, nr0_sq, stage1_dir)

E_vmc_minus_eb = E_vmc .- eb.(h_vmc) ./ 2
E_dmc_minus_eb = E_dmc .- eb.(h_dmc) ./ 2

pgfplotsx()
p = scatter(h_vmc, E_vmc_minus_eb; yerror=err_vmc, label="VMC",
            xlabel=L"h/r_0", ylabel=L"E/N - \varepsilon_b/2 \, [\hbar^2/(mr_0^2)]",
            marker=:circle, markersize=8, color=:dodgerblue, framestyle=:box, legend=:topright)
scatter!(p, h_dmc, E_dmc_minus_eb; yerror=err_dmc, label="DMC", marker=:diamond, markersize=8, color=:darkorange)
isfinite(E_single) && hline!(p, [E_single]; color=:black, linestyle=:dash,
                              label=L"nr_0^2=%$(nr0_sq/2)\ \mathrm{(single\ layer)}")

savefig(p, joinpath(stage_dir, "data", "plots", "Fig1_combined_N$(num_part)_nr0sq$(nr0_sq).pdf"))
display(p)
println("Saved combined Fig1 plot")