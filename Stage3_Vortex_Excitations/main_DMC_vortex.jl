# main_DMC_bilayer.jl
using DelimitedFiles, Plots, LaTeXStrings, Statistics, ProgressMeter, Base.Threads

include(normpath(joinpath(@__DIR__, "src", "utils.jl")))
include(normpath(joinpath(@__DIR__, "src", "jastrow.jl")))
include(normpath(joinpath(@__DIR__, "src", "energy.jl")))
include(normpath(joinpath(@__DIR__, "src", "metropolis.jl")))
include(normpath(joinpath(@__DIR__, "src", "dmc.jl")))
include(normpath(joinpath(@__DIR__, "src", "shooting_method.jl")))
include(normpath(joinpath(@__DIR__, "src", "observables.jl")))
include(normpath(joinpath(@__DIR__, "scripts", "dmc_extrapolation.jl")))

"""
    load_vmc_density(vmc_data_path, N, h, d, lA, lB)

Loads the VMC-sampled n(x,y) grids for both layers and the shared xy_bins
produced by run_vmc_vortex_sweep.jl, for the mixed-estimator extrapolation
against this DMC run's density at the same (N, h, d, lA, lB).
"""
function load_vmc_density(vmc_data_path::String, N::Int, h::Float64, d::Float64,
                           lA::Float64, lB::Float64)
    d_str     = round(d, digits=3)
    nA_path   = joinpath(vmc_data_path, "nxyA_N$(N)_h$(h)_d$(d_str)_lA$(lA)_lB$(lB).txt")
    nB_path   = joinpath(vmc_data_path, "nxyB_N$(N)_h$(h)_d$(d_str)_lA$(lA)_lB$(lB).txt")
    bins_path = joinpath(vmc_data_path, "xy_bins_N$(N).txt")

    isfile(nA_path)   || error("VMC density file not found: $nA_path")
    isfile(nB_path)   || error("VMC density file not found: $nB_path")
    isfile(bins_path) || error("VMC bins file not found: $bins_path")

    n_xy_A_vmc  = readdlm(nA_path, '\t', Float64)
    n_xy_B_vmc  = readdlm(nB_path, '\t', Float64)
    xy_bins_vmc = vec(readdlm(bins_path, '\t', Float64))

    return n_xy_A_vmc, n_xy_B_vmc, xy_bins_vmc
end

# ── Parameters ────────────────────────────────────────────────────────────────
N             = 60
nr0sq         = 1.0
L             = sqrt(N / nr0sq)
num_walkers   = 200
quadratic     = true           # quadratic vs linear short-time propagator
lA            = 1.0            # vortex circulation for layer A
lB            = 1.0            # vortex circulation for layer B
prop_tag      = quadratic ? "quad" : "lin"   # goes into every filename (see note below)
total_time    = 10.0
equil_time    = 2.0
τ_factor      = 0.1

num_bins      = 100            # density/current grid resolution for dmc()'s n_xy_A/n_xy_B
n_r_bins      = 30             # radial bins for the core-centered profiles
min_count_warn = 20            # grid points below which a radial bin is not trustworthy

h_vals_to_run = if !isempty(ARGS)
    [parse(Float64, ARGS[1])]
else
    collect(range(0.3, 1.1, step = 0.7))
end

r_min     = 1e-6
Δ_shoot   = 1e-4
tol_shoot = 1e-10

# ── R_match (Stage I, AA/BB — h-independent) ────────────────────────────────
vmc_path_stage1 = joinpath(@__DIR__, "..", "Stage1_2D_Dipol_System", "data", "results", "VMC",
           "vmc_N$(N ÷ 2)_nr0sq$(nr0sq / 2).txt")
vmc_data_stage1 = readdlm(vmc_path_stage1, '\t', String, skipstart=1)
R_match = parse(Float64, vmc_data_stage1[1, 4])

# d-sweep, identical convention to run_vmc_vortex_sweep.jl so VMC and DMC
# results land on the same d grid and are directly comparable.
# d_vals = collect(0.0:0.05*(L/4):L/4)
d_vals = [L/4]
n_d    = length(d_vals)

# ── Output directories ──────────────────────────────────────────────────────
plots_path       = joinpath(@__DIR__, "data", "results", "DMC", "plots")
vmc_density_path = joinpath(@__DIR__, "..", "Stage3_Vortex_Excitations", "data", "results", "VortexSweep", "data")
mixed_data_path  = joinpath(@__DIR__, "data", "results", "DMC", "mixed_density")
extr_data_path   = joinpath(@__DIR__, "data", "results", "DMC", "extrapolated_density")
radial_data_path = joinpath(@__DIR__, "data", "results", "DMC", "radial")
gr_data_path     = joinpath(@__DIR__, "data", "results", "DMC", "pair_correlations")
for p in (plots_path, mixed_data_path, extr_data_path, radial_data_path, gr_data_path)
    mkpath(p)
end

"""
    radial_density_profile(n_xy, xy_bins, x0, y0, L; n_r_bins=30, θ0=nothing, Δθ=π/4)

Core-centered radial density profile.

Returns `(r_eff, n_r, counts)`, where `r_eff[k]` is the mean radius of the grid
points that fell into bin `k`, NOT the geometric bin center. In a shell the
area-weighted mean radius exceeds the midpoint, and the offset is largest in the
innermost bins — exactly where the n(r) ~ r^2 core law is fitted. Plotting
against bin centers therefore flattens the log-log slope artificially.

If `θ0` is given, only grid points whose polar angle about (x0,y0) lies within
`Δθ` of `θ0` are accumulated. Used to separate the direction pointing at the
other layer's core from the complementary direction, so that the induced
modulation is not washed out by the angular average.

Duplicated from run_vmc_vortex_sweep.jl — apply any change to BOTH copies, or
the VMC and DMC profiles stop being comparable point by point. Worth moving
into observables.jl once both scripts stabilize.
"""
function radial_density_profile(n_xy::Matrix{Float64}, xy_bins::Vector{Float64},
                                 x0::Float64, y0::Float64, L::Float64;
                                 n_r_bins::Int=30,
                                 θ0::Union{Nothing,Float64}=nothing,
                                 Δθ::Float64=π/4)
    num_bins_grid = length(xy_bins)
    dr            = (L/2) / n_r_bins
    r_sum         = zeros(Float64, n_r_bins)
    r_rad         = zeros(Float64, n_r_bins)
    r_count       = zeros(Int,     n_r_bins)

    @inbounds for i in 1:num_bins_grid, j in 1:num_bins_grid
        dx = get_periodic_difference(xy_bins[i], x0, L)
        dy = get_periodic_difference(xy_bins[j], y0, L)
        r  = sqrt(dx^2 + dy^2)
        r < L/2 || continue

        if θ0 !== nothing
            dθ = mod(atan(dy, dx) - θ0 + π, 2π) - π
            abs(dθ) <= Δθ || continue
        end

        bin_index = clamp(Int(floor(r / dr)) + 1, 1, n_r_bins)
        r_sum[bin_index]   += n_xy[i, j]
        r_rad[bin_index]   += r
        r_count[bin_index] += 1
    end

    r_eff = [r_count[k] > 0 ? r_rad[k] / r_count[k] : NaN for k in 1:n_r_bins]
    n_r   = [r_count[k] > 0 ? r_sum[k] / r_count[k] : NaN for k in 1:n_r_bins]
    return r_eff, n_r, Float64.(r_count)
end

# ── Main loop over h ─────────────────────────────────────────────────────────
for h in h_vals_to_run

    println("\n" * "="^100)
    println("BILAYER DMC PRODUCTION RUN   h = $h r₀,  lA = $lA,  lB = $lB,  N = $N,  nr0sq = $nr0sq,  Δτ_factor = $τ_factor,  prop = $prop_tag")
    println("="^100)

    # ── R0_opt + E_ref_initial (Stage II VMC results, production subfolder) ──
    vmc_path_stage2 = joinpath(@__DIR__, "..", "Stage2_Bilayer_Dipolar_Bosons", "data", "results", "VMC", "production",
               "VMC_Stage2_N$(N)_nr0sq$(nr0sq)_h$(h).txt")
    vmc_data_stage2 = readdlm(vmc_path_stage2, '\t', String, skipstart=1)
    R0_opt        = parse(Float64, vmc_data_stage2[1, 4])
    E_ref_initial = parse(Float64, vmc_data_stage2[1, 5]) * N

    Constants = calculate_constants(L, R_match)

    # ── Δτ_ref from the move_all delta sweep (find_delta_bilayer.jl output) ─
    Δτ_ref = 10^(-3)

    # ── Build fAB (cubic-spline interlayer Jastrow derivatives) ─────────────
    _, _, _, itp_up, itp_upp, _ =
        build_fAB(h, R0_opt, r_min, Δ_shoot, tol_shoot, nr0sq, N)

    N_half = N ÷ 2

    println("  R0_opt        = $R0_opt")
    println("  R_match       = $R_match")
    println("  E_ref_initial = $(E_ref_initial/N)")
    println("  Δτ_ref        = $Δτ_ref")

    Δτ        = Δτ_ref * τ_factor
    num_steps = max(10^4, round(Int, total_time / Δτ))
    num_equil = max(2000, round(Int, equil_time / Δτ))

    results_path = joinpath(@__DIR__, "data", "results", "DMC",
                   "DMC_Stage3_N$(N)_nr0sq$(nr0sq)_h$(h)_lA$(lA)_lB$(lB)_$(prop_tag).txt")
    mkpath(dirname(results_path))

    # Per-h accumulators over the d sweep, for the E(d) summary plot
    d_recorded  = Float64[]
    E_per_N     = Float64[]
    E_per_N_err = Float64[]

    open(results_path, "w") do io
        println(io, "d\tE_dmc\tE_dmc_err\tnum_steps\tnum_equil\tΔτ")

        # ── d loop — sequential (see thread-nesting note above) ─────────────
        for idx in 1:n_d
            d     = d_vals[idx]
            d_str = round(d, digits=3)
            suffix = "N$(N)_h$(h)_d$(d_str)_lA$(lA)_lB$(lB)_$(prop_tag)"

            x_vortex_A, y_vortex_A = L/2, L/2
            x_vortex_B, y_vortex_B = L/2 + d, L/2

            println("\n═══ h=$h | d=$d | lA=$lA | lB=$lB | Δτ=$(round(Δτ, sigdigits=3)) | " *
                    "steps=$num_steps | walkers=$num_walkers ═══")

            # ── Initial configuration: random placement ─────────────────────
            x0, y0 = random_initial_config(N, L, "Uniform")
            xA_init, yA_init = x0[1:N_half],     y0[1:N_half]
            xB_init, yB_init = x0[N_half+1:end], y0[N_half+1:end]

            E_dmc, E_dmc_err, E_history,
            n_xy_A, n_xy_B, xy_bins,
            gAA_r, gBB_r, gAB_r, g_total_r, r_vals =
                dmc(
                    xA_init, yA_init, xB_init, yB_init,
                    num_walkers, N, num_steps,
                    Δτ, L, h,
                    lA, lB, R_match, Constants,
                    x_vortex_A, y_vortex_A, x_vortex_B, y_vortex_B,
                    R0_opt, itp_up, itp_upp,
                    E_ref_initial, num_walkers;
                    num_equil   = num_equil,
                    quadratic   = quadratic,
                    plot_energy = false,
                    num_bins    = num_bins
                )

            # Block averaging
            block_sizes = [10, 20, 30, 40, 50, 100, 150, 200,
                           300, 400, 500, 600, 700, 800, 900, 1000,
                           1100, 1200, 1300, 1400, 1500, 1600, 1700,
                           1800, 1900, 2000]
            sigmas  = [blocking_statistics(E_history, B)[2] for B in block_sizes]
            plateau = detect_plateau(block_sizes, sigmas, window_size=4, rtol=0.02)
            avg_E, σ = blocking_statistics(E_history, plateau)

            println("E/N = $(round(avg_E/N, digits=5)) ± $(round(σ/N, digits=5))")

            println(io, "$(d)\t$(avg_E)\t$(σ)\t$(num_steps)\t$(num_equil)\t$(Δτ)")

            push!(d_recorded, d)
            push!(E_per_N, avg_E/N)
            push!(E_per_N_err, σ/N)

            # ── Persist the raw mixed estimator and the DMC grid ─────────────
            # The mixed density is the object the extrapolations are built from;
            # without it on disk nothing downstream can be redone without a rerun.
            writedlm(joinpath(mixed_data_path, "nxyA_mixed_$(suffix).txt"), n_xy_A)
            writedlm(joinpath(mixed_data_path, "nxyB_mixed_$(suffix).txt"), n_xy_B)
            writedlm(joinpath(mixed_data_path, "xy_bins_$(suffix).txt"),   xy_bins)

            # ── Pair correlations (returned by dmc(), previously discarded) ──
            writedlm(joinpath(gr_data_path, "g_r_$(suffix).txt"),
                     hcat(r_vals, gAA_r, gBB_r, gAB_r, g_total_r))

            p_gr = plot(r_vals, gAA_r; label = L"g_{AA}(r)", xlabel = L"r", ylabel = L"g(r)",
                        title = L"h=%$h,\ d=%$(d_str),\ \ell_A=%$lA,\ \ell_B=%$lB",
                        framestyle = :box, grid = true, gridalpha = 0.25, legend = :bottomright)
            plot!(p_gr, r_vals, gBB_r;    label = L"g_{BB}(r)")
            plot!(p_gr, r_vals, gAB_r;    label = L"g_{AB}(r)")
            plot!(p_gr, r_vals, g_total_r; label = L"g_{tot}(r)", linestyle = :dash)
            savefig(p_gr, joinpath(plots_path, "g_r_$(suffix).pdf"))

            # ── n(x,y) heatmap, mixed estimator ──────────────────────────────
            p_density = plot(
                heatmap(xy_bins, xy_bins, n_xy_A'; title="Layer A", c=:viridis, aspect_ratio=:equal),
                heatmap(xy_bins, xy_bins, n_xy_B'; title="Layer B", c=:viridis, aspect_ratio=:equal),
                layout = (1, 2), size = (900, 400),
                plot_title = L"n(x,y)\ \mathrm{(mixed)},\ \ h=%$h,\ d=%$(d_str),\ \ell_A=%$lA,\ \ell_B=%$lB"
            )
            scatter!(p_density[1], [x_vortex_A], [y_vortex_A]; marker=:xcross, ms=8, mc=:red, label="core A")
            scatter!(p_density[2], [x_vortex_B], [y_vortex_B]; marker=:xcross, ms=8, mc=:red, label="core B")
            savefig(p_density, joinpath(plots_path, "density_DMC_$(suffix).pdf"))

            # ── VMC densities, for the extrapolation and the overlay ────────
            n_xy_A_vmc, n_xy_B_vmc, xy_bins_vmc = load_vmc_density(vmc_density_path, N, h, d, lA, lB)

            # Float equality on the grids is fragile; compare with a tolerance.
            if length(xy_bins_vmc) != length(xy_bins) || !all(isapprox.(xy_bins_vmc, xy_bins; atol=1e-10))
                @warn "h=$h, d=$d: VMC and DMC xy_bins differ — check num_bins/L consistency " *
                      "between run_vmc_vortex_sweep.jl and this DMC run."
            end

            # ── Core-centered radial profiles: mixed vs VMC ─────────────────
            r_A,     n_r_A,     counts_A     = radial_density_profile(n_xy_A,     xy_bins, x_vortex_A, y_vortex_A, L; n_r_bins=n_r_bins)
            r_B,     n_r_B,     counts_B     = radial_density_profile(n_xy_B,     xy_bins, x_vortex_B, y_vortex_B, L; n_r_bins=n_r_bins)
            r_A_vmc, n_r_A_vmc, _            = radial_density_profile(n_xy_A_vmc, xy_bins, x_vortex_A, y_vortex_A, L; n_r_bins=n_r_bins)
            r_B_vmc, n_r_B_vmc, _            = radial_density_profile(n_xy_B_vmc, xy_bins, x_vortex_B, y_vortex_B, L; n_r_bins=n_r_bins)

            writedlm(joinpath(radial_data_path, "radial_$(suffix).txt"),
                     hcat(r_A, n_r_A, n_r_A_vmc, counts_A, r_B, n_r_B, n_r_B_vmc, counts_B))

            # The sign of (mixed - VMC) near r ~ d says whether the true ground
            # state has the induced modulation more or less strongly than the
            # ansatz, i.e. which way the extrapolation will push it.
            p_radial = plot(
                r_A, n_r_A;
                seriestype = :scatter, label = "layer A (mixed)",
                xlabel = L"r", ylabel = L"n(r)",
                title  = L"h=%$h,\ d=%$(d_str),\ \ell_A=%$lA,\ \ell_B=%$lB",
                xscale = :log10, yscale = :log10,
                legend = :bottomright, framestyle = :box,
                grid = true, gridalpha = 0.25
            )
            scatter!(p_radial, r_B,     n_r_B;     label = "layer B (mixed)")
            plot!(   p_radial, r_A_vmc, n_r_A_vmc; label = "layer A (VMC)", linestyle = :dash, seriestype = :line)
            plot!(   p_radial, r_B_vmc, n_r_B_vmc; label = "layer B (VMC)", linestyle = :dot,  seriestype = :line)
            vline!(  p_radial, [d]; label = L"r = d", linestyle = :dashdot, linecolor = :gray)
            savefig(p_radial, joinpath(plots_path, "radial_density_DMC_$(suffix).pdf"))

            if any(counts_A[1:3] .< min_count_warn) || any(counts_B[1:3] .< min_count_warn)
                @warn "h=$h, d=$d: fewer than $min_count_warn grid points in one of the innermost radial bins " *
                      "(counts_A[1:3]=$(counts_A[1:3]), counts_B[1:3]=$(counts_B[1:3])). " *
                      "Increase num_bins or decrease n_r_bins before reading a core exponent off these points."
            end

            # ── Angular decomposition: toward vs away from the other core ────
            if d > 0
                θ_AB = atan(get_periodic_difference(y_vortex_B, y_vortex_A, L),
                            get_periodic_difference(x_vortex_B, x_vortex_A, L))

                r_A_tw, n_A_tw, c_A_tw = radial_density_profile(n_xy_A, xy_bins, x_vortex_A, y_vortex_A, L; n_r_bins=n_r_bins, θ0=θ_AB)
                r_A_aw, n_A_aw, c_A_aw = radial_density_profile(n_xy_A, xy_bins, x_vortex_A, y_vortex_A, L; n_r_bins=n_r_bins, θ0=θ_AB + π)
                r_B_tw, n_B_tw, c_B_tw = radial_density_profile(n_xy_B, xy_bins, x_vortex_B, y_vortex_B, L; n_r_bins=n_r_bins, θ0=θ_AB + π)
                r_B_aw, n_B_aw, c_B_aw = radial_density_profile(n_xy_B, xy_bins, x_vortex_B, y_vortex_B, L; n_r_bins=n_r_bins, θ0=θ_AB)

                writedlm(joinpath(radial_data_path, "radial_wedge_$(suffix).txt"),
                         hcat(r_A_tw, n_A_tw, c_A_tw, r_A_aw, n_A_aw, c_A_aw,
                              r_B_tw, n_B_tw, c_B_tw, r_B_aw, n_B_aw, c_B_aw))

                p_wedge = plot(
                    r_A_tw, n_A_tw;
                    seriestype = :scatter, label = "layer A, toward core B",
                    xlabel = L"r", ylabel = L"n(r)",
                    title  = L"\mathrm{angular\ sectors},\ h=%$h,\ d=%$(d_str)",
                    legend = :bottomright, framestyle = :box,
                    grid = true, gridalpha = 0.25
                )
                scatter!(p_wedge, r_A_aw, n_A_aw; label = "layer A, away")
                scatter!(p_wedge, r_B_tw, n_B_tw; label = "layer B, toward core A")
                scatter!(p_wedge, r_B_aw, n_B_aw; label = "layer B, away")
                vline!(  p_wedge, [d]; label = L"r = d", linestyle = :dashdot, linecolor = :gray)
                savefig(p_wedge, joinpath(plots_path, "radial_wedge_$(suffix).pdf"))
            end

            # ── DMC/VMC mixed-estimator extrapolation of n(x,y), Eq. B.14 ──
            n_xy_A_extr1, n_xy_A_extr2 = extrapolate_mixed(n_xy_A, n_xy_A_vmc)
            n_xy_B_extr1, n_xy_B_extr2 = extrapolate_mixed(n_xy_B, n_xy_B_vmc)

            writedlm(joinpath(extr_data_path, "nxyA_extr1_$(suffix).txt"), n_xy_A_extr1)
            writedlm(joinpath(extr_data_path, "nxyA_extr2_$(suffix).txt"), n_xy_A_extr2)
            writedlm(joinpath(extr_data_path, "nxyB_extr1_$(suffix).txt"), n_xy_B_extr1)
            writedlm(joinpath(extr_data_path, "nxyB_extr2_$(suffix).txt"), n_xy_B_extr2)

            # Diagnostics for the two pathologies seen at small h: the linear
            # form can go negative where n_VMC is small, and the quadratic form
            # is a pointwise ratio with no normalization constraint.
            cell_area = (xy_bins[2] - xy_bins[1])^2
            for (tag, m) in (("A_extr1", n_xy_A_extr1), ("A_extr2", n_xy_A_extr2),
                             ("B_extr1", n_xy_B_extr1), ("B_extr2", n_xy_B_extr2))
                n_neg = count(<(0), m)
                norm  = sum(m) * cell_area
                println("    $(tag): negative cells = $n_neg,  ∫n dA = $(round(norm, digits=4)) (expected $(N ÷ 2))")
            end

            p_extr = plot(
                heatmap(xy_bins, xy_bins, n_xy_A_extr1'; title="Layer A, linear extr.",    c=:viridis, aspect_ratio=:equal),
                heatmap(xy_bins, xy_bins, n_xy_B_extr1'; title="Layer B, linear extr.",    c=:viridis, aspect_ratio=:equal),
                heatmap(xy_bins, xy_bins, n_xy_A_extr2'; title="Layer A, quadratic extr.", c=:viridis, aspect_ratio=:equal),
                heatmap(xy_bins, xy_bins, n_xy_B_extr2'; title="Layer B, quadratic extr.", c=:viridis, aspect_ratio=:equal),
                layout = (2, 2), size = (900, 800),
                plot_title = L"n_{extr}(x,y),\ \ h=%$h,\ d=%$(d_str),\ \ell_A=%$lA,\ \ell_B=%$lB"
            )
            savefig(p_extr, joinpath(plots_path, "density_extrapolated_$(suffix).pdf"))
        end
    end

    # ── E_dmc(d)/N summary plot for this h ─────────────────────────────────
    p_E = plot(d_recorded ./ (L/2), E_per_N;
               yerror     = E_per_N_err,
               marker     = :circle,
               linewidth  = 1.5,
               label      = "DMC",
               xlabel     = L"d\ /\ (L/2)",
               ylabel     = L"E_{DMC}/N\ \ (\varepsilon_0)",
               title      = L"h = %$h,\ \ell_A = %$lA,\ \ell_B = %$lB",
               legend     = :topright,
               framestyle = :box,
               grid       = true,
               gridalpha  = 0.25)
    savefig(p_E, joinpath(plots_path, "dmc_energy_vs_offset_N$(N)_h$(h)_lA$(lA)_lB$(lB)_$(prop_tag).pdf"))

    println("Saved: $results_path")
    println("Plots saved to $plots_path")
end

println("\nAll bilayer DMC production runs completed.")