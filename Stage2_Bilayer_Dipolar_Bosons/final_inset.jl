# Stage2_Bilayer_Dipolar_Bosons/final_inset.jl
#
# Standalone — does not need main_VMC.jl/main_DMC.jl to have run in the
# same session. Requires both VMC and DMC pipelines already saved for
# every h in h_vals. Inset figure only: E/N (+ tail) and ε_b/2 vs h/r₀,
# VMC and DMC together.

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

h_range = collect(LinRange(minimum(h_vals), maximum(h_vals), 200))

pgfplotsx()
p = scatter(h_vmc, E_vmc; yerror=err_vmc, label="VMC",
            xlabel=L"h/r_0", ylabel=L"E/N \, [\hbar^2/(mr_0^2)]",
            marker=:circle, markersize=8, color=:dodgerblue, framestyle=:box, legend=:topleft)
scatter!(p, h_dmc, E_dmc; yerror=err_dmc, label="DMC", marker=:diamond, markersize=8, color=:darkorange)
plot!(p, h_range, eb.(h_range) ./ 2; label=L"\varepsilon_b/2", linewidth=2, color=:black)
hline!(p, [0.0]; color=:black, linestyle=:dash, label="")
isfinite(E_single) && hline!(p, [E_single]; color=:red, linestyle=:dash,
                              label=L"nr_0^2=%$(nr0_sq/2)\ \mathrm{(single\ layer)}")

savefig(p, joinpath(stage_dir, "data", "plots", "Inset_combined_N$(num_part)_nr0sq$(nr0_sq).pdf"))
display(p)
println("Saved combined inset plot")