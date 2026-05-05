# scripts/plot_results.jl
# Expects from main.jl: num_part, nr0_sq
# Reads results from data/results/ and generates plots

using Plots, LaTeXStrings, PGFPlotsX

num_part = 30
nr0_sq = 16.0

L = sqrt(num_part / nr0_sq)
nr0_sq_32 = num_part * nr0_sq^(3/2)

println("Reading results for N = $num_part, nr0^2 = $nr0_sq")

pgfplotsx()

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
        line == "R_match\tEnergy\tError_Energy\tDrift_Energy\tError_Drift\tLaplacian_Energy\tError_Laplacian" && continue
        line == "R_opt\tE_opt"    && continue
        vals = parse.(Float64, split(line, "\t"))
        if occursin("Coarse", section)
            push!(R_coarse, vals[1]); push!(E_coarse, vals[2]);push!(error_coarse, vals[3])
            push!(R_coarse_drift, vals[1]);push!(E_coarse_drift, vals[4]);push!(error_coarse_drift, vals[5])
            push!(R_coarse_laplacian, vals[1]); push!(E_coarse_laplacian, vals[6]);push!(error_coarse_laplacian, vals[7])
        elseif occursin("Fine", section)
            push!(R_fine, vals[1]); push!(E_fine, vals[2]);push!(error_fine, vals[3])
            push!(R_fine_drift, vals[1]);push!(E_fine_drift, vals[4]);push!(error_fine_drift, vals[5])
            push!(R_fine_laplacian, vals[1]); push!(E_fine_laplacian, vals[6]);push!(error_fine_laplacian, vals[7])
        elseif occursin("Optimal", section)
            R_opt[] = vals[1]
            E_opt[] = vals[2]
        end
    end
end

println("Loaded $(length(R_coarse)) coarse and $(length(R_fine)) fine sweep points")
println("Optimal R_match from file: ", R_opt[])

# ------------------------------------------
# Plot: R_match sweep
# ------------------------------------------
default(
    fontfamily = "Computer Modern",
    titlefontsize = 14,
    guidefontsize = 12,
    tickfontsize  = 10,
    legendfontsize = 9,
    grid = true,
    gridalpha = 0.3,
    framestyle = :box,
    dpi = 600
)

col_std = "#3182bd"
col_drift = "#31a354"
col_lap = "#d95f02"

p1 = plot(
    xlabel  = L"R_{\mathrm{match}}",
    ylabel  = L"(E/N) \cdot (nr_0^2)^{-3/2}",
    title   = L"R_{\mathrm{match}}\ \mathrm{Optimization},\ nr_0^2 = %$(Int(nr0_sq)),\ N = %$(num_part)",
    legend  = :topright,
    dpi = 600
)

# Coarse sweep — transparent, no legend clutter
plot!(R_coarse, E_coarse ./ nr0_sq_32,
      yerror = error_coarse ./ nr0_sq_32,
      label = "", color = col_std, alpha = 0.3,
      marker = :circle, markersize = 3, markerstrokewidth = 0)

plot!(R_coarse_drift, E_coarse_drift ./ nr0_sq_32,
      yerror = error_coarse_drift ./ nr0_sq_32,
      label = "", color = col_drift, alpha = 0.3,
      marker = :diamond, markersize = 3, markerstrokewidth = 0)

plot!(R_coarse_laplacian, E_coarse_laplacian ./ nr0_sq_32,
      yerror = error_coarse_laplacian ./ nr0_sq_32,
      label = "", color = col_lap, alpha = 0.3,
      marker = :square, markersize = 3, markerstrokewidth = 0)

# Fine sweep — solid, with legend
plot!(R_fine, E_fine ./ nr0_sq_32,
      yerror = error_fine ./ nr0_sq_32,
      label = L"E_{\mathrm{loc}}", color = col_std, linewidth = 2,
      marker = :circle, markersize = 5, markerstrokewidth = 0)

plot!(R_fine_drift, E_fine_drift ./ nr0_sq_32,
      yerror = error_fine_drift ./ nr0_sq_32,
      label = L"E_{\mathrm{drift}}", color = col_drift, linewidth = 2,
      marker = :diamond, markersize = 5, markerstrokewidth = 0)

plot!(R_fine_laplacian, E_fine_laplacian ./ nr0_sq_32,
      yerror = error_fine_laplacian ./ nr0_sq_32,
      label = L"E_{\mathrm{lap}}", color = col_lap, linewidth = 2,
      marker = :square, markersize = 5, markerstrokewidth = 0)

# Optimal R_match
vline!([R_opt[]],
       linestyle = :dash, linewidth = 1.5, color = :black, alpha = 0.8,
       label = L"R_{\mathrm{opt}} = %$(round(R_opt[], digits=4))")

savefig(joinpath(@__DIR__, "..", "data", "plots",
        "plot_Rmatch_N$(num_part)_nr0sq$(nr0_sq).pdf"))

println("Saved plot to data/plots/plot_Rmatch_N$(num_part)_nr0sq$(nr0_sq).pdf")
display(p1)
