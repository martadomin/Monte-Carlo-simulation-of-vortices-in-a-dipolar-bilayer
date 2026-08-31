# scripts/run_vmc_vortex_sweep.jl
using DelimitedFiles, Plots, LaTeXStrings, Statistics, ProgressMeter, Base.Threads, Printf

include(normpath(joinpath(@__DIR__, "src", "utils.jl")))
include(normpath(joinpath(@__DIR__, "src", "jastrow.jl")))
include(normpath(joinpath(@__DIR__, "src", "energy.jl")))
include(normpath(joinpath(@__DIR__, "src", "metropolis.jl")))
include(normpath(joinpath(@__DIR__, "src", "observables.jl")))
include(normpath(joinpath(@__DIR__, "src", "shooting_method.jl")))

# ── Parameters ────────────────────────────────────────────────────────────
N        = 60
num_part = N
nr0sq    = 1.0
L        = sqrt(N / nr0sq)
l        = 1.0                              # vortex circulation
h_vals   = range(0.7, 0.7, step = 0.1)       # sweep over interlayer separations
num_steps_production = 10^7

r_min     = 1e-6
Δ_shoot   = 1e-4
tol_shoot = 1e-10

# R_match doesn't depend on h — load once, outside the h loop
nr0sq_half_str = @sprintf("%.1f", nr0sq/2)
rmatch_path_stage1 = joinpath(@__DIR__, "..", "Stage1_2D_Dipol_System", "data", "results", "VMC",
                       "vmc_N$(N÷2)_nr0sq$(nr0sq_half_str).txt")
rmatch_vals = readdlm(rmatch_path_stage1, '\t', String, skipstart=1)
R_match = parse(Float64, rmatch_vals[1, 4])
Constants = calculate_constants(L, R_match)

# Output directory for the plots (defined once, outside the h loop)
plots_path = joinpath(@__DIR__, "data", "results", "VortexSweep", "plots")
mkpath(plots_path)

for h in h_vals

    # R0 and the AB splines DO depend on h — reload/rebuild per h
    nr0sq_str = @sprintf("%.1f", nr0sq)
    r0_path_stage2 = joinpath(@__DIR__, "..", "Stage2_Bilayer_Dipolar_Bosons", "data", "results", "VMC", "production",
               "VMC_Stage2_N$(N)_nr0sq$(nr0sq_str)_h$(h).txt")
    r0_vals = readdlm(r0_path_stage2, '\t', String, skipstart=1)
    R0     = parse(Float64, r0_vals[1, 4])
    E0     = parse(Float64, r0_vals[1, 5])
    err_E0 = parse(Float64, r0_vals[1, 6])
    println("Radius:", R0)
    println("Bilayer Energy result VMC:", E0, " ± ", err_E0)

    _, _, itp_u, itp_up, itp_upp, _ =
        build_fAB(h, R0, r_min, Δ_shoot, tol_shoot, nr0sq, N)


    d_vals = collect(0.0:0.05*(L/4):L/4)

    # ── Baseline is the Stage II scalar E0 ± err_E0 (single number, no
    # d-dependence) — at l=0 the Stage III ansatz reduces exactly to the
    # plain Stage II wavefunction (f_0(r) ≡ 1), so E0 is the same physical
    # quantity for ANY l, and is subtracted identically here regardless of
    # what l this run uses. No special-casing l==0 is needed anymore.

    results_path = joinpath(@__DIR__, "data", "results", "VortexSweep",
                             "vortex_energy_vs_offset_N$(num_part)_h$(h)_l$(l).txt")

    mkpath(dirname(results_path))

    n_d = length(d_vals)

    # Accumulators for the post-loop plot — fresh for each h, filled in d-order after threading
    d_recorded       = Float64[]
    E_per_N          = Float64[]   # E_vortex/N, baseline-subtracted against Stage II E0
    E_per_N_err      = Float64[]
    E_drift          = Float64[]   # E_vortex_drift/N, baseline-subtracted against Stage II E0
    E_drift_err      = Float64[]
    E_laplacian      = Float64[]   # E_vortex_laplacian/N, baseline-subtracted against Stage II E0
    E_laplacian_err  = Float64[]

    # ── Pre-allocated per-d result slots (written by independent threads) ──
    E_tot_arr           = fill(NaN, n_d)
    E_err_arr           = fill(NaN, n_d)
    E_tot_drift_arr     = fill(NaN, n_d)
    E_err_drift_arr     = fill(NaN, n_d)
    E_tot_laplacian_arr = fill(NaN, n_d)
    E_err_laplacian_arr = fill(NaN, n_d)
    E_kin_arr           = fill(NaN, n_d)
    E_int_arr           = fill(NaN, n_d)
    acc_arr             = fill(NaN, n_d)

    # Each d is an independent VMC chain — no warm-start chaining across d
    # once threaded (there's no well-defined "previous d" when chains run
    # concurrently), so every thread starts from its own random configuration.
    Threads.@threads for idx in 1:n_d
        d = d_vals[idx]
        x_vortex_A, y_vortex_A = L/2, L/2
        x_vortex_B, y_vortex_B = L/2 + d, L/2

        x0, y0 = random_initial_config(num_part, L, "Uniform")
        x_A_i, y_A_i = x0[1:num_part÷2], y0[1:num_part÷2]
        x_B_i, y_B_i = x0[num_part÷2+1:end], y0[num_part÷2+1:end]

        delta_tuned, x_A_t, y_A_t, x_B_t, y_B_t = tune_delta(
            x_A_i, y_A_i, x_B_i, y_B_i,
            L, R_match, Constants,
            x_vortex_A, y_vortex_A, x_vortex_B, y_vortex_B, l,
            R0, itp_u
        )

        energies, energies_drift, energies_laplacian,
        E_tot, E_sq, E_tot_drift, E_tot_laplacian,
        E_kin, E_int, acc_ratio, x_coord, y_coord =
            metropolis(num_part, num_steps_production, delta_tuned, L, h,
                       R_match, R0, itp_u, itp_up, itp_upp, Constants;
                       x_vortex_A=x_vortex_A, y_vortex_A=y_vortex_A,
                       x_vortex_B=x_vortex_B, y_vortex_B=y_vortex_B, l=l,
                       x_A_init=x_A_t, y_A_init=y_A_t,
                       x_B_init=x_B_t, y_B_init=y_B_t,
                       progress=false)

        # Block-average all three estimators independently, matching run_vmc.jl:
        # each gets its own plateau, since drift/Laplacian have different
        # variance/autocorrelation behavior than the local (std) estimator.
        block_sizes = [10, 20, 30, 40, 50, 100, 150, 200,
                       300, 400, 500, 600, 700, 800, 900, 1000,
                       1100, 1200, 1300, 1400, 1500, 1600, 1700,
                       1800, 1900, 2000, 3000, 5000, 10000, 20000, 40000]

        sigmas           = [blocking_statistics(energies, B)[2]            for B in block_sizes]
        sigmas_drift     = [blocking_statistics(energies_drift, B)[2]      for B in block_sizes]
        sigmas_laplacian = [blocking_statistics(energies_laplacian, B)[2]  for B in block_sizes]

        plateau_std       = detect_plateau(block_sizes, sigmas,           window_size=4, rtol=0.05)
        plateau_drift     = detect_plateau(block_sizes, sigmas_drift,     window_size=4, rtol=0.05)
        plateau_laplacian = detect_plateau(block_sizes, sigmas_laplacian, window_size=4, rtol=0.05)

        E_tot_b,           E_err_b            = blocking_statistics(energies,           plateau_std)
        E_tot_drift_b,     E_err_drift_b      = blocking_statistics(energies_drift,     plateau_drift)
        E_tot_laplacian_b, E_err_laplacian_b  = blocking_statistics(energies_laplacian, plateau_laplacian)

        E_tot_arr[idx]           = E_tot_b
        E_err_arr[idx]           = E_err_b
        E_tot_drift_arr[idx]     = E_tot_drift_b
        E_err_drift_arr[idx]     = E_err_drift_b
        E_tot_laplacian_arr[idx] = E_tot_laplacian_b
        E_err_laplacian_arr[idx] = E_err_laplacian_b
        E_kin_arr[idx]           = E_kin
        E_int_arr[idx]           = E_int
        acc_arr[idx]             = acc_ratio

        println("[thread $(Threads.threadid())] h=$h, d=$d done: " *
                "E_std=$(E_tot_b/num_part) ± $(E_err_b/num_part) (block=$plateau_std)  acc=$acc_ratio")
    end

    # ── Baseline subtraction against the Stage II scalar E0 ± err_E0 ───────
    # Applied identically to all three estimators (std, drift, laplacian), so
    # that the Eq. 70 three-estimator cross-check is performed on genuinely
    # comparable (baseline-subtracted) quantities, not a mix of subtracted
    # and raw values. E0 is d-independent, so it is subtracted the same way
    # from every d point, for any l (no l==0 special-casing needed since E0
    # comes from an external Stage II file, not from this script's own output).
    E_vortex_arr               = E_tot_arr           .- E0
    E_vortex_err_arr           = sqrt.(E_err_arr.^2           .+ err_E0^2)
    E_vortex_drift_arr         = E_tot_drift_arr     .- E0
    E_vortex_drift_err_arr     = sqrt.(E_err_drift_arr.^2     .+ err_E0^2)
    E_vortex_laplacian_arr     = E_tot_laplacian_arr .- E0
    E_vortex_laplacian_err_arr = sqrt.(E_err_laplacian_arr.^2 .+ err_E0^2)

    # ── Sequential I/O and plot accumulation, in d order ────────────────────
    open(results_path, "w") do io
        println(io, "d\tE_tot\tE_tot_err\tE_vortex\tE_vortex_err\t" *
                     "E_drift\tE_drift_err\tE_laplacian\tE_laplacian_err\tE_kin\tE_int\tacceptance")

        for idx in 1:n_d
            d = d_vals[idx]
            println(io, "$(d)\t$(E_tot_arr[idx])\t$(E_err_arr[idx])\t" *
                        "$(E_vortex_arr[idx])\t$(E_vortex_err_arr[idx])\t" *
                        "$(E_vortex_drift_arr[idx])\t$(E_vortex_drift_err_arr[idx])\t" *
                        "$(E_vortex_laplacian_arr[idx])\t$(E_vortex_laplacian_err_arr[idx])\t" *
                        "$(E_kin_arr[idx])\t$(E_int_arr[idx])\t$(acc_arr[idx])")

            push!(d_recorded, d / (L/2))
            push!(E_per_N, E_vortex_arr[idx] / num_part)
            push!(E_per_N_err, E_vortex_err_arr[idx] / num_part)
            push!(E_drift, E_vortex_drift_arr[idx] / num_part)
            push!(E_drift_err, E_vortex_drift_err_arr[idx] / num_part)
            push!(E_laplacian, E_vortex_laplacian_arr[idx] / num_part)
            push!(E_laplacian_err, E_vortex_laplacian_err_arr[idx] / num_part)
        end
    end

    # ── Plot E_vortex(d)/N with error bars, all three estimators ──────────
    p = plot(d_recorded, E_per_N;
              yerror     = E_per_N_err,
              marker     = :circle,
              linewidth  = 1.5,
              label      = "local (std)",
              xlabel     = L"d\ /\ (L/2)",
              ylabel     = L"E_{vortex}/N\ \ (\varepsilon_0)",
              title      = L"h = %$h,\ \ell = %$l",
              legend     = :topright,
              framestyle = :box,
              grid       = true,
              gridalpha  = 0.25)

    scatter!(p, d_recorded, E_drift;     yerror = E_drift_err,     label = "drift")
    scatter!(p, d_recorded, E_laplacian; yerror = E_laplacian_err, label = "laplacian")

    display(p)
    fig_path = joinpath(plots_path, "vortex_energy_vs_offset_N$(num_part)_h$(h)_l$(l).pdf")
    savefig(p, fig_path)

    println("\nVortex offset sweep complete for h=$h. Results saved to $results_path")
    println("Plot saved to $fig_path")
end