# Stage3_Vortex_Excitations/scripts/final_plot_vortex_DMC.jl
#
# E_vortex(d)/N vs d — DMC only. Same shape as final_plot_vortex_VMC.jl.

using Plots, LaTeXStrings

d_found, E_vals, err_vals = Float64[], Float64[], Float64[]
for d in d_vals
    path = result_path(stage_dir, "DMC_vortex", (N=num_part, nr0sq=nr0_sq, h=h, lA=lA, lB=lB, d=d))
    if !isfile(path)
        @warn "No saved DMC_vortex result for d=$d — skipping"
        continue
    end
    r = load_run(path; run=1).result
    push!(d_found, d / (L/2)); push!(E_vals, r.E_extrapolated / num_part); push!(err_vals, r.E_extrapolated_err / num_part)
end

p = plot(d_found, E_vals; yerror=err_vals, marker=:diamond, linewidth=1.5, color=:darkorange, label="DMC",
         xlabel=L"d\ /\ (L/2)", ylabel=L"E_{vortex}/N\ (\varepsilon_0)",
         title=L"h = %$h,\ \ell_A=%$lA,\ \ell_B=%$lB", legend=:topright, framestyle=:box)

display(p)
savefig(p, joinpath(stage_dir, "data", "plots", "vortex_energy_vs_d_DMC_N$(num_part)_nr0sq$(nr0_sq)_h$(h).pdf"))
println("Saved DMC vortex offset plot")