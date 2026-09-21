# Stage2_Bilayer_Dipolar_Bosons/scripts/final_inset_DMC.jl
#
# Inset figure: E/N (+ tail) and ε_b/2 vs h/r₀ — DMC only.
# Expects: N, nr0sq, h_vals, stage_dir, stage1_dir

using Plots, LaTeXStrings

h_found, E_corr_vals, err_vals = load_h_sweep_dmc(N, nr0sq, h_vals, stage_dir)
eb = exact_binding_energy_interp(stage_dir)
E_single = single_layer_reference(N, nr0sq, stage1_dir)

h_range = collect(LinRange(minimum(h_found), maximum(h_found), 200))

pgfplotsx()
p = scatter(h_found, E_corr_vals; yerror=err_vals, label="DMC",
            xlabel=L"h/r_0", ylabel=L"E/N \, [\hbar^2/(mr_0^2)]",
            marker=:diamond, markersize=8, color=:darkorange, framestyle=:box, legend=:topleft)
plot!(p, h_range, eb.(h_range) ./ 2; label=L"\varepsilon_b/2", linewidth=2, color=:black)
hline!(p, [0.0]; color=:black, linestyle=:dash, label="")
isfinite(E_single) && hline!(p, [E_single]; color=:red, linestyle=:dash,
                              label=L"nr_0^2=%$(nr0sq/2)\ \mathrm{(single\ layer)}")

mkpath(joinpath(stage_dir, "data", "plots"))
savefig(p, joinpath(stage_dir, "data", "plots", "Inset_DMC_N$(N)_nr0sq$(nr0sq).pdf"))
display(p)