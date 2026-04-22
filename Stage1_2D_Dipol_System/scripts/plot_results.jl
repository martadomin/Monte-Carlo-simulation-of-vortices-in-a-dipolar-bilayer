# scripts/plot_results.jl
# Expects from main.jl: num_part, nr0_sq
# Reads results from data/results/ and generates plots

using Plots, LaTeXStrings
gr()

# ------------------------------------------
# Load results
# ------------------------------------------
rmatch_path = joinpath(@__DIR__, "..", "data", "results",
              "Rmatch_sweep_N$(num_part)_nr0sq$(nr0_sq).txt")

vmc_path = joinpath(@__DIR__, "..", "data", "results",
           "vmc_N$(num_part)_nr0sq$(nr0_sq).txt")

# Parse R_match sweep file
R_coarse, E_coarse = Float64[], Float64[]
R_fine,   E_fine   = Float64[], Float64[]
R_opt, E_opt       = 0.0, 0.0

open(rmatch_path, "r") do io
    section = ""
    for line in eachline(io)
        startswith(line, "#") && (section = line; continue)
        isempty(strip(line))  && continue
        line == "R_match\tEnergy" && continue
        line == "R_opt\tE_opt"    && continue
        vals = parse.(Float64, split(line, "\t"))
        if occursin("Coarse", section)
            push!(R_coarse, vals[1]); push!(E_coarse, vals[2])
        elseif occursin("Fine", section)
            push!(R_fine, vals[1]);   push!(E_fine, vals[2])
        elseif occursin("Optimal", section)
            R_opt = vals[1];          E_opt = vals[2]
        end
    end
end

# ------------------------------------------
# Plot 1: R_match sweep
# ------------------------------------------
p1 = plot(R_coarse, E_coarse,
          label=L"Coarse sweep",
          xlabel=L"R_{\mathrm{match}}",
          ylabel=L"E/N \cdot (nr_0^2)^{-3/2}",
          title="R_match optimization, nr0^2 = $(nr0_sq), N = $(num_part)",
          marker=:circle,
          linewidth=2)
plot!(R_fine, E_fine,
      label=L"Fine sweep",
      marker=:circle,
      linewidth=2)
vline!([R_opt], linestyle=:dash, color=:red,
       label="R_opt = $(round(R_opt, digits=4))")

savefig(joinpath(@__DIR__, "..", "data", "results",
        "plot_Rmatch_N$(num_part)_nr0sq$(nr0_sq).pdf"))

# ------------------------------------------
# Plot 2: VMC energy vs paper
# ------------------------------------------
# Paper fit: E/N = a1*(nr0^2)^(3/2) + a2*(nr0^2)^(5/4) + a3*(nr0^2)^(1/2)
a1, a2, a3 = 4.536, 4.38, 1.2  # gas phase coefficients from Astrakharchik 2007
nr0_range  = collect(LinRange(1.0, 300.0, 500))
E_paper    = @. a1*nr0_range^(3/2) + a2*nr0_range^(5/4) + a3*nr0_range^(1/2)
E_paper_normalized = E_paper ./ nr0_range.^(3/2)

# Read VMC result
vmc_data = readdlm(vmc_path, '\t', Float64, skipstart=1)
E_vmc    = vmc_data[1, 5]  # E/N/(nr0^2)^(3/2) column

p2 = plot(nr0_range, E_paper_normalized,
          label=L"DMC fit (Astrakharchik 2007)",
          xlabel=L"nr_0^2",
          ylabel=L"E/N \cdot (nr_0^2)^{-3/2}",
          title="VMC vs DMC, N = $(num_part)",
          linewidth=2,
          color=:blue)
          
scatter!([nr0_sq], [E_vmc],
         label="VMC, N = $(num_part)",
         marker=:circle,
         markersize=8,
         color=:red)

savefig(joinpath(@__DIR__, "..", "data", "results",
        "plot_vmc_vs_paper_N$(num_part)_nr0sq$(nr0_sq).pdf"))

display(p1)
display(p2)