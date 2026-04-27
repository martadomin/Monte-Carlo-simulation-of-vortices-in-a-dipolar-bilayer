# scripts/plot_results.jl
# Expects from main.jl: num_part, nr0_sq
# Reads results from data/results/ and generates plots

using Plots, LaTeXStrings

num_part = 30
nr0_sq = 16.0
nr0_sq_32 = num_part * nr0_sq^(3/2)

println("Reading results for N = $num_part, nr0^2 = $nr0_sq")

gr()

# ------------------------------------------
# Load results
# ------------------------------------------
rmatch_path = joinpath(@__DIR__, "..", "data", "sweep_results",
              "Rmatch_sweep_N$(num_part)_nr0sq$(nr0_sq).txt")


# Parse R_match sweep file
R_coarse, E_coarse, error_coarse = Float64[], Float64[], Float64[]
R_coarse_drift, E_coarse_drift, error_coarse_drift = Float64[], Float64[], Float64[]
R_coarse_laplacian, E_coarse_laplacian, error_coarse_laplacian = Float64[], Float64[], Float64[]
R_fine, E_fine, error_fine = Float64[], Float64[], Float64[]
R_fine_drift, E_fine_drift, error_fine_drift = Float64[], Float64[], Float64[]
R_fine_laplacian, E_fine_laplacian, error_fine_laplacian = Float64[], Float64[], Float64[]
R_opt, E_opt = Ref(0.0), Ref(0.0)

open(rmatch_path, "r") do io
    section = ""
    for line in eachline(io)
        startswith(line, "#") && (section = line; continue)
        isempty(strip(line))  && continue
        line == "R_match\tEnergy\tDrift_Energy\tLaplacian_Energy" && continue
        line == "R_opt\tE_opt"    && continue
        vals = parse.(Float64, split(line, "\t"))
        if occursin("Coarse", section)
            push!(R_coarse, vals[1]); push!(E_coarse, vals[2])
            push!(R_coarse_drift, vals[1]);push!(E_coarse_drift, vals[3])
            push!(R_coarse_laplacian, vals[1]); push!(E_coarse_laplacian, vals[4])
        elseif occursin("Fine", section)
            push!(R_fine, vals[1]); push!(E_fine, vals[2])
            push!(R_fine_drift, vals[1]);push!(E_fine_drift, vals[3])
            push!(R_fine_laplacian, vals[1]); push!(E_fine_laplacian, vals[4])
        elseif occursin("Optimal", section)
            R_opt[] = vals[1]
            E_opt[] = vals[2]
        end
    end
end

println("Loaded $(length(R_coarse)) coarse and $(length(R_fine)) fine sweep points")
println("Optimal R_match from file: ", R_opt[])

# ------------------------------------------
# Plot 1: R_match sweep
# ------------------------------------------
# Add error bars to the plot
p1 = plot(R_coarse, E_coarse ./ nr0_sq_32,
          label=L"\textrm{Coarse \; sweep}",
          xlabel=L"R_{\mathrm{match}}",
          ylabel=L"E/N \cdot (nr_0^2)^{-3/2}",
          title=L"R_{\mathrm{match}}\textrm{\; optimization, \;} nr_0^2 = %$(nr0_sq),\ N = %$(num_part)",
          marker=:circle,
          linewidth=2,
          legend=:topleft)

plot!(R_fine, E_fine ./ nr0_sq_32,  
      label=L"\textrm{Fine \; sweep}",
      marker=:circle,
      linewidth=2)
      
plot!(R_coarse_drift, E_coarse_drift ./ nr0_sq_32,
      label=L"\textrm{Coarse \; sweep (drift)}",
      marker=:diamond,
      linewidth=2)
    
plot!(R_fine_drift, E_fine_drift ./ nr0_sq_32,
      label=L"\textrm{Fine \; sweep (drift)}",
      marker=:diamond,
      linewidth=2)

plot!(R_coarse_laplacian, E_coarse_laplacian ./ nr0_sq_32,
      label=L"\textrm{Coarse \; sweep (Laplacian)}",
      marker=:square,
      linewidth=2)

plot!(R_fine_laplacian, E_fine_laplacian ./ nr0_sq_32,
      label=L"\textrm{Fine \; sweep (Laplacian)}",
      marker=:square,
      linewidth=2)

vline!([R_opt[]], linestyle=:dash, color=:red,
       label=L"R_{\mathrm{opt}} = %$(round(R_opt[], digits=4))")


savefig(joinpath(@__DIR__, "..", "data", "plots",
        "plot_Rmatch_N$(num_part)_nr0sq$(nr0_sq).pdf"))
println("Saved R_match sweep plot to: ", joinpath("data", "plots", "plot_Rmatch_N$(num_part)_nr0sq$(nr0_sq).pdf"))

display(p1)
