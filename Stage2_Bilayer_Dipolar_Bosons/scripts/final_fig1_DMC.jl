# Stage2_Bilayer_Dipolar_Bosons/scripts/final_fig1_DMC.jl
#
# Fig1: E/N - ε_b/2 vs h/r₀ — DMC only.
# Expects: N, nr0sq, h_vals, stage_dir, stage1_dir

using Plots, LaTeXStrings

h_found, E_corr_vals, err_vals = load_h_sweep_dmc(N, nr0sq, h_vals, stage_dir)
eb = exact_binding_energy_interp(stage_dir)
E_single = single_layer_reference(N, nr0sq, stage1_dir)

E_minus_eb_vals = E_corr_vals .- eb.(h_found) ./ 2

pgfplotsx()
p = scatter(h_found, E_minus_eb_vals; yerror=err_vals, label="DMC",
            xlabel=L"h/r_0", ylabel=L"E/N - \varepsilon_b/2 \, [\hbar^2/(mr_0^2)]",
            marker=:diamond, markersize=8, color=:darkorange, framestyle=:box, legend=:topright)
isfinite(E_single) && hline!(p, [E_single]; color=:black, linestyle=:dash,
                              label=L"nr_0^2=%$(nr0sq/2)\ \mathrm{(single\ layer)}")

savefig(p, joinpath(stage_dir, "data", "plots", "Fig1_DMC_N$(N)_nr0sq$(nr0sq).pdf"))
display(p)