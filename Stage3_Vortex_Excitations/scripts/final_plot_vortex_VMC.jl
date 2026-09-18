# Stage3_Vortex_Excitations/scripts/final_plot_vortex_VMC.jl
#
# E_vortex(d)/N vs d — VMC only. Expects: num_part, nr0_sq, h, lA, lB,
# d_vals, L, stage_dir

using Plots, LaTeXStrings

d_found, E_vals, err_vals = Float64[], Float64[], Float64[]
for d in d_vals
    path = result_path(stage_dir, "VMC_vortex", (N=num_part, nr0sq=nr0_sq, h=h, lA=lA, lB=lB, d=d))
    if !isfile(path)
        @warn "No saved VMC_vortex result for d=$d — skipping"
        continue
    end
    r = load_run(path; run=1).result
    push!(d_found, d / (L/2)); push!(E_vals, r.avg_energy / num_part); push!(err_vals, r.sigma / num_part)
end

p = plot(d_found, E_vals; yerror=err_vals, marker=:circle, linewidth=1.5, label="VMC",
         xlabel=L"d\ /\ (L/2)", ylabel=L"E_{vortex}/N\ (\varepsilon_0)",
         title=L"h = %$h,\ \ell_A=%$lA,\ \ell_B=%$lB", legend=:topright, framestyle=:box)

display(p)
savefig(p, joinpath(stage_dir, "data", "plots", "vortex_energy_vs_d_VMC_N$(num_part)_nr0sq$(nr0_sq)_h$(h).pdf"))
println("Saved VMC vortex offset plot")