# scripts/plot_R0_sweep.jl
#
# Plots the R0 optimization sweep (coarse + fine) for a single h value,
# Stage II bilayer — local, drift, and Laplacian kinetic-energy estimators.
# Style analogous to Stage I's plot_results.jl.
#
# Expects globals from main_VMC_bilayer.jl: N, nr0sq, h, L
# Reads: data/sweep_results/R0_sweep_Stage2_N$(N)_nr0sq$(nr0sq)_h$(h).txt
# File columns: R0  E/N  Error  eps_b  E_drift/N  Error_drift  E_lap/N  Error_lap

println("Reading R0 sweep for N = $N, nr0² = $nr0sq, h = $h")

pgfplotsx()

# ------------------------------------------
# Load results
# ------------------------------------------
sweep_path = joinpath(@__DIR__, "..", "data", "sweep_results",
                      "R0_sweep_Stage2_N$(N)_nr0sq$(nr0sq)_h$(h).txt")

R_coarse, E_coarse, error_coarse, eb_coarse = Float64[], Float64[], Float64[], Float64[]
E_coarse_drift, error_coarse_drift = Float64[], Float64[]
E_coarse_laplacian, error_coarse_laplacian = Float64[], Float64[]
R_fine, E_fine, error_fine, eb_fine = Float64[], Float64[], Float64[], Float64[]
E_fine_drift, error_fine_drift = Float64[], Float64[]
E_fine_laplacian, error_fine_laplacian = Float64[], Float64[]
R0_opt_file = err_opt_file = E_opt_file = energy_b_file = NaN

open(sweep_path, "r") do io
    section = ""
    for line in eachline(io)
        s = line

        # ── Comment lines ────────────────────────────────────────
        if startswith(s, "#")
            # Only update section on section-marker lines, not column headers
            if occursin("Coarse sweep", s) || occursin("Fine sweep", s) || occursin("Optimal", s)
                section = s
            end
            continue
        end

        isempty(strip(s)) && continue

        # ── Optimal block: key = value lines (no tabs) ───────────
        if occursin("=", s) && !occursin("\t", s)
            k, v = strip.(split(s, "=", limit=2))
            k == "R0_opt"   && (R0_opt_file   = parse(Float64, v))
            k == "E_opt"    && (E_opt_file    = parse(Float64, v))
            k == "err_opt"  && (err_opt_file  = parse(Float64, v))
            k == "energy_b" && (energy_b_file = parse(Float64, v))
            continue
        end

        # ── Data lines: tab-separated ────────────────────────────
        vals = try
            parse.(Float64, split(s, "\t"))
        catch
            continue
        end

        if occursin("Coarse", section)
            push!(R_coarse,     vals[1]); push!(E_coarse,      vals[2])
            push!(error_coarse, vals[3]); push!(eb_coarse,     vals[4])
            if length(vals) >= 8
                push!(E_coarse_drift,      vals[5]); push!(error_coarse_drift,      vals[6])
                push!(E_coarse_laplacian,  vals[7]); push!(error_coarse_laplacian,  vals[8])
            end
        elseif occursin("Fine", section)
            push!(R_fine,     vals[1]); push!(E_fine,      vals[2])
            push!(error_fine, vals[3]); push!(eb_fine,     vals[4])
            if length(vals) >= 8
                push!(E_fine_drift,      vals[5]); push!(error_fine_drift,      vals[6])
                push!(E_fine_laplacian,  vals[7]); push!(error_fine_laplacian,  vals[8])
            end
        end
    end
end

println("Loaded $(length(R_coarse)) coarse and $(length(R_fine)) fine sweep points")
println("Optimal R0 from file: ", R0_opt_file, "   E/N = ", E_opt_file, " ± ", err_opt_file)

isempty(R_coarse) && isempty(R_fine) &&
    error("No sweep points parsed from $sweep_path — check file format/section markers.")

have_drift_lap = length(E_fine_drift) == length(R_fine) && !isempty(E_fine_drift)

# ------------------------------------------
# Plot: R0 sweep — local, drift, laplacian
# ------------------------------------------
default(
    fontfamily     = "Computer Modern",
    titlefontsize  = 14,
    guidefontsize  = 12,
    tickfontsize   = 10,
    legendfontsize = 9,
    grid           = true,
    gridalpha      = 0.3,
    framestyle     = :box,
    dpi            = 600
)

col_std   = "#3182bd"
col_drift = "#31a354"
col_lap   = "#d95f02"

p1 = plot(
    xlabel = L"R_0 \, [r_0]",
    ylabel = L"E/N \, [\hbar^2/(mr_0^2)]",
    title  = L"R_0\ \mathrm{Optimization},\ nr_0^2 = %$(nr0sq),\ N = %$(N),\ h = %$(round(h, digits=3))",
    legend = :topright,
    dpi    = 600
)

# Coarse sweep — transparent, no legend clutter
plot!(p1, R_coarse, E_coarse,
      yerror = error_coarse,
      label = "", color = col_std, alpha = 0.3,
      marker = :circle, markersize = 3, markerstrokewidth = 0)

if have_drift_lap
    plot!(p1, R_coarse, E_coarse_drift,
          yerror = error_coarse_drift,
          label = "", color = col_drift, alpha = 0.3,
          marker = :diamond, markersize = 3, markerstrokewidth = 0)

    plot!(p1, R_coarse, E_coarse_laplacian,
          yerror = error_coarse_laplacian,
          label = "", color = col_lap, alpha = 0.3,
          marker = :square, markersize = 3, markerstrokewidth = 0)
end

# Fine sweep — solid, with legend
plot!(p1, R_fine, E_fine,
      yerror = error_fine,
      label = L"E_{\mathrm{loc}}", color = col_std, linewidth = 2,
      marker = :circle, markersize = 5, markerstrokewidth = 0)

if have_drift_lap
    plot!(p1, R_fine, E_fine_drift,
          yerror = error_fine_drift,
          label = L"E_{\mathrm{drift}}", color = col_drift, linewidth = 2,
          marker = :diamond, markersize = 5, markerstrokewidth = 0)

    plot!(p1, R_fine, E_fine_laplacian,
          yerror = error_fine_laplacian,
          label = L"E_{\mathrm{lap}}", color = col_lap, linewidth = 2,
          marker = :square, markersize = 5, markerstrokewidth = 0)
end

# Optimal R0
if !isnan(R0_opt_file)
    vline!(p1, [R0_opt_file],
           linestyle = :dash, linewidth = 1.5, color = :black, alpha = 0.8,
           label = L"R_{0,\mathrm{opt}} = %$(round(R0_opt_file, digits=4))")
end

plots_dir = joinpath(@__DIR__, "..", "data", "plots")
mkpath(plots_dir)
savefig(p1, joinpath(plots_dir, "plot_R0_sweep_N$(N)_nr0sq$(nr0sq)_h$(h).pdf"))

println("Saved plot to data/plots/plot_R0_sweep_N$(N)_nr0sq$(nr0sq)_h$(h).pdf")
display(p1)