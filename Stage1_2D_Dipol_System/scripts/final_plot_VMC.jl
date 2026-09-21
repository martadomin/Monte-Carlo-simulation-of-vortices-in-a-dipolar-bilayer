# Stage1_2D_Dipol_System/scripts/final_plot_VMC.jl
#
# Loads every VMC result saved by run_vmc.jl for this num_part, plots
# E/N vs nr0_sq against the Astrakharchik 2007 fit. Expects: num_part,
# nr0_sq_vals, stage_dir.

using Plots, LaTeXStrings

plot_convergence_analysis = false   # set true to also show the blocking-σ(B) curves

E_vmc, error_E_vmc = Float64[], Float64[]
for nr0_sq in nr0_sq_vals
    p_result = result_path(stage_dir, "VMC", (N=num_part, nr0sq=nr0_sq))
    loaded_run = load_run(p_result; run=1)
    push!(E_vmc, loaded_run.result.avg_energy / (num_part * nr0_sq^(3/2)))
    push!(error_E_vmc, loaded_run.result.sigma / (num_part * nr0_sq^(3/2)))
end

a1, a2, a3 = 4.536, 4.38, 1.2
E_tail = tail_energy(num_part)
nr0_range = collect(LinRange(10.0, 300.0, 500))
E_paper_normalized = (@. a1*nr0_range^(3/2) + a2*nr0_range^(5/4) + a3*nr0_range^(1/2)) ./ nr0_range.^(3/2)

gr()
p = plot(nr0_range, E_paper_normalized,
          label=L"\mathrm{DMC\ fit\ (Astrakharchik\ 2007)}",
          xlabel=L"nr_0^2", ylabel=L"E/N \cdot (nr_0^2)^{-3/2}",
          title=L"\mathrm{VMC},\ N = %$(num_part)",
          linewidth=1, color=:black, framestyle=:box, legend=:topright)
scatter!(nr0_sq_vals, E_vmc .+ E_tail, yerror=error_E_vmc,
         label=L"\mathrm{VMC,\ }N = %$(num_part)", marker=:circle, markersize=6, color=:red)

mkpath(joinpath(stage_dir, "data", "plots"))
savefig(joinpath(stage_dir, "data", "plots", "plot_vmc_N$(num_part).pdf"))
println("Saved VMC plot")
display(p)

if plot_convergence_analysis
    for nr0_sq in nr0_sq_vals
        run = load_run(result_path(stage_dir, "VMC", (N=num_part, nr0sq=nr0_sq)); run=1)
        r = run.result
        p_conv = plot_convergence(r.block_sizes, r.sigmas, r.sigmas_drift, r.sigmas_laplacian,
                                   r.plateau_std, r.plateau_drift, r.plateau_laplacian;
                                   title="Convergence, nr0²=$nr0_sq")
        mkpath(joinpath(stage_dir, "data", "plots"))
        display(p_conv)
        savefig(p_conv, joinpath(stage_dir, "data", "plots", "convergence_N$(num_part)_nr0sq$(nr0_sq).pdf"))
    end
end