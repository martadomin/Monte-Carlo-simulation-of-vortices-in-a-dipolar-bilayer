using DelimitedFiles, Plots, LaTeXStrings

num_part = 30
nr0_sq_sim = [16.0, 32.0, 48.0, 64.0, 96.0, 128.0, 196.0, 256.0, 384.0, 512.0, 768.0, 1024.0]
E_vmc = Float64[]
error_E_vmc = Float64[]

for nr0_sq in nr0_sq_sim
        E_tail = @. π * sqrt(nr0_sq) / sqrt(num_part)
        vmc_path = joinpath(@__DIR__, "..", "data", "results",
                "vmc_N$(num_part)_nr0sq$(nr0_sq).txt")
        vmc_data = readdlm(vmc_path, '\t', Float64, skipstart=1)
        nr0_sq_32 = num_part * nr0_sq^(3/2)
        push!(E_vmc, vmc_data[1, 5] / nr0_sq_32 + E_tail)        
        push!(error_E_vmc, vmc_data[1, 6] / nr0_sq_32) 
end

# Paper fit: E/N = a1*(nr0^2)^(3/2) + a2*(nr0^2)^(5/4) + a3*(nr0^2)^(1/2)
a1, a2, a3 = 4.536, 4.38, 1.2
nr0_range  = collect(LinRange(1.0, 1050.0, 500))
E_paper    = @. a1*nr0_range^(3/2) + a2*nr0_range^(5/4) + a3*nr0_range^(1/2)
E_paper_normalized = E_paper ./ nr0_range.^(3/2)

pgfplotsx()

p2 = plot(nr0_range, E_paper_normalized,
          label=L"DMC\ \mathrm{fit\ (Astrakharchik\ 2007)}",
          xlabel=L"nr_0^2",
          ylabel=L"E/N \cdot (nr_0^2)^{-3/2}",
          title=L"VMC\ \mathrm{vs\ DMC},\ N = %$(num_part)",
          linewidth=2,
          color=:blue,
          framestyle=:box,
          legend=:topright)

scatter!(nr0_sq_sim, E_vmc,
         yerror=error_E_vmc,
         label=L"VMC,\ N = %$(num_part)",
         marker=:circle,
         markersize=6,
         color=:red)

savefig(joinpath(@__DIR__, "..", "data", "plots",
        "plot_vmc_vs_paper_N$(num_part).pdf"))
println("Saved final plot")
display(p2)