nr0_sq_sim = [16.0, 32.0, 48.0, 64.0, 96.0, 128.0, 196.0, 256.0, 384.0, 512.0, 768.0, 1024.0]
E_vmc = Float64[]
error_E_vmc = Float64[]

for nr0_sq in nr0_sq_sim
    L = sqrt(num_part / nr0_sq)
    vmc_path = joinpath(@__DIR__, "..", "data", "results",
            "vmc_N$(num_part)_nr0sq$(nr0_sq).txt")

    # Read VMC result
    vmc_data = readdlm(vmc_path, '\t', Float64, skipstart=1)
    push!(E_vmc, vmc_data[1, 5]/num_part)  # Energy per particle
    push!(error_E_vmc, vmc_data[1, 6]/num_part)
end

# ------------------------------------------
# Plot 2: VMC energy vs paper
# ------------------------------------------
# Paper fit: E/N = a1*(nr0^2)^(3/2) + a2*(nr0^2)^(5/4) + a3*(nr0^2)^(1/2)

a1, a2, a3 = 4.536, 4.38, 1.2  # gas phase coefficients from Astrakharchik 2007
nr0_range  = collect(LinRange(1.0, 300.0, 500))
E_paper    = @. a1*nr0_range^(3/2) + a2*nr0_range^(5/4) + a3*nr0_range^(1/2)
E_paper_normalized = E_paper ./ nr0_range.^(3/2)


p2 = plot(nr0_range, E_paper_normalized,
          label=L"DMC fit (Astrakharchik 2007)",
          xlabel=L"nr_0^2",
          ylabel=L"E/N \cdot (nr_0^2)^{-3/2}",
          title="VMC vs DMC, N = $(num_part)",
          linewidth=2,
          color=:blue)
          
scatter!(nr0_sq_sim, E_vmc ./ nr0_sq_sim.^(3/2),
         label="VMC, N = $(num_part)",
         marker=:circle,
         markersize=8,
         color=:red)

savefig(joinpath(@__DIR__, "..", "data", "results",
        "plot_vmc_vs_paper_N$(num_part)_nr0sq$(nr0_sq).pdf"))

display(p2)