# Stage3_Vortex_Excitations/final_plot_vortex.jl
#
# Standalone — combined VMC+DMC, E_vortex(d)/N vs d. Requires both
# main_VMC.jl and main_DMC.jl already run for this sweep.

include(joinpath(@__DIR__, "..", "common", "src", "CommonCore.jl"))
include(joinpath(@__DIR__, "src", "wavefunction.jl"))
include(joinpath(@__DIR__, "config.jl"))
using Plots, LaTeXStrings

d_vmc, E_vmc, err_vmc = Float64[], Float64[], Float64[]
for d in d_vals
    path = result_path(stage_dir, "VMC_vortex", (N=num_part, nr0sq=nr0_sq, h=h, lA=lA, lB=lB, d=d))
    isfile(path) || continue
    r = load_run(path; run=1).result
    push!(d_vmc, d/(L/2)); push!(E_vmc, r.avg_energy/num_part); push!(err_vmc, r.sigma/num_part)
end

d_dmc, E_dmc, err_dmc = Float64[], Float64[], Float64[]
for d in d_vals
    path = result_path(stage_dir, "DMC_vortex", (N=num_part, nr0sq=nr0_sq, h=h, lA=lA, lB=lB, d=d))
    isfile(path) || continue
    r = load_run(path; run=1).result
    push!(d_dmc, d/(L/2)); push!(E_dmc, r.E_extrapolated/num_part); push!(err_dmc, r.E_extrapolated_err/num_part)
end

p = scatter(d_vmc, E_vmc; yerror=err_vmc, marker=:circle, label="VMC", color=:dodgerblue,
            xlabel=L"d\ /\ (L/2)", ylabel=L"E_{vortex}/N\ (\varepsilon_0)",
            title=L"h = %$h,\ \ell_A=%$lA,\ \ell_B=%$lB", legend=:topright, framestyle=:box)
scatter!(d_dmc, E_dmc; yerror=err_dmc, marker=:diamond, label="DMC", color=:darkorange)

display(p)
savefig(p, joinpath(stage_dir, "data", "plots", "vortex_energy_vs_d_combined_N$(num_part)_nr0sq$(nr0_sq)_h$(h).pdf"))
println("Saved combined vortex offset plot")