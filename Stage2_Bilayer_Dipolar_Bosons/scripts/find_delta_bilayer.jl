"""
find_delta_bilayer.jl

Standalone script to find the optimal delta for move_all VMC in the bilayer system.
Reads R_match from the Stage I VMC results (AA/BB) and R0_opt + E_ref from
the Stage II VMC results (AB), using the same readdlm convention as run_dmc.jl.
Does NOT call optimize_R0.jl or run_VMC_bilayer.jl.

Sets on exit:
  delta_opt_moveall, Δτ_DMC
"""

using Random, DelimitedFiles, Plots, LaTeXStrings, ProgressMeter, Dierckx, Base.Threads

include(normpath(joinpath(@__DIR__, "..", "src", "utils.jl")))
include(normpath(joinpath(@__DIR__, "..", "src", "jastrow.jl")))
include(normpath(joinpath(@__DIR__, "..", "src", "energy.jl")))
include(normpath(joinpath(@__DIR__, "..", "src", "metropolis.jl")))
include(normpath(joinpath(@__DIR__, "..", "src", "shooting_method.jl")))
include(normpath(joinpath(@__DIR__, "..", "src", "observables.jl")))

# ── Parameters ─────────────────────────────────────────────────────────────────
N      = 60
nr0sq  = 1.0
h_vals      = range(0.3, 1.5, step = 0.1)
L      = sqrt(N / nr0sq)

r_min     = 1e-6
Δ_shoot   = 1e-4
tol_shoot = 1e-10

num_sweep_steps = 10^5
delta_vals = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.12, 0.15, 0.20, 0.25]

@threads for h in h_vals
    # ── Read R_match (Stage I, AA/BB) and R0_opt + E_ref (Stage II, AB) ────────────
    # Mirrors the convention used in scripts/run_dmc.jl

    vmc_path_stage1 = joinpath(@__DIR__, "..", "..", "Stage1_2D_Dipol_System", "data", "results", "VMC",
            "vmc_N$(N ÷ 2)_nr0sq$(nr0sq / 2).txt")
    vmc_data_stage1 = readdlm(vmc_path_stage1, '\t', String, skipstart=1)
    R_match = parse(Float64, vmc_data_stage1[1, 4])

    vmc_path_stage2 = joinpath(@__DIR__, "..", "data", "results", "VMC", "production",
            "VMC_Stage2_N$(N)_nr0sq$(nr0sq)_h$(h).txt")
    vmc_data_stage2 = readdlm(vmc_path_stage2, '\t', String, skipstart=1)
    R0_opt        = parse(Float64, vmc_data_stage2[1, 4])

    println("="^70)
    println("DELTA SWEEP (move_all)  h = $h r₀")
    println("="^70)
    println("  R0_opt        = $R0_opt")
    println("  R_match       = $R_match")
    println("  L             = $L\n")

    Constants = calculate_constants(L, R_match)

    # ── Build fAB ─────────────────────────────────────────────────────────────────
    _, _, itp_u, itp_up, itp_upp, _ =
        build_fAB(h, R0_opt, r_min, Δ_shoot, tol_shoot, nr0sq, N)

    # ── Initial configuration ──────────────────────────────────────────────────────
    x_coord, y_coord = random_initial_config(N, L, "Uniform")
    x_A = x_coord[1:N÷2];     y_A = y_coord[1:N÷2]
    x_B = x_coord[N÷2+1:end]; y_B = y_coord[N÷2+1:end]

    # ── Sweep ──────────────────────────────────────────────────────────────────────
    acceptance_vals = Float64[]

    for δ in delta_vals
        _, _, _, _, _, _, _, _, _, acc, _, _ = metropolis(
            N, num_sweep_steps, δ, L, h, R_match, R0_opt,
            itp_u, itp_up, itp_upp, Constants;
            x_A_init = x_A, y_A_init = y_A,
            x_B_init = x_B, y_B_init = y_B,
            progress = true,
            move_all = true
        )
        push!(acceptance_vals, acc)
        println("  δ = $(rpad(string(δ), 6)) → acceptance = $(round(acc * 100, digits=2))%")
    end

    # ── Find optimal delta via interpolation (acceptance = 50%) ────────────────────
    # Interpolate δ as a function of acceptance (sorted by acceptance, ascending),
    # then evaluate at acceptance = 0.5. Interpolation is in log10(δ) since the
    # sweep is log-spaced and acceptance vs log10(δ) is closer to linear.

    sort_idx   = sortperm(acceptance_vals)
    acc_sorted = acceptance_vals[sort_idx]
    logδ_sorted = log10.(delta_vals[sort_idx])

    @assert minimum(acc_sorted) <= 0.5 <= maximum(acc_sorted) "Acceptance never crosses 50% within the sampled δ range"

    # Dierckx's Spline1D requires strictly increasing knots and default k=3 (cubic)
    # needs at least 4 points; falls back to a lower degree automatically otherwise.
    itp_logδ = Spline1D(acc_sorted, logδ_sorted; k = min(3, length(acc_sorted) - 1))
    delta_opt_moveall = 10.0^itp_logδ(0.5)
    acc_opt           = 0.5

    Δτ_DMC = delta_opt_moveall^2   # D = 1/2, Δτ = δ²/(2D) = δ²

    println("\n  Optimal δ (interpolated) = $delta_opt_moveall  (acceptance = $(round(acc_opt * 100, digits=2))%)")
    println("  Δτ_DMC                   = $Δτ_DMC")

    # ── Plot ───────────────────────────────────────────────────────────────────────
    p = plot(delta_vals, acceptance_vals .* 100;
            xlabel  = L"\delta",
            ylabel  = "Acceptance (%)",
            title   = "move_all delta sweep  N=$N, h=$h r₀",
            marker  = :circle,
            xscale  = :log10,
            label   = "acceptance",
            legend  = :topright)
    hline!(p, [50.0]; linestyle = :dash, color = :red,   label = "50%")
    vline!(p, [delta_opt_moveall]; linestyle = :dot, color = :black, label = "optimal δ")
    display(p)
    readline()

    # ── Save ───────────────────────────────────────────────────────────────────────
    out_path = joinpath(@__DIR__, "..", "data", "sweep_results",
            "delta_moveall_N$(N)_nr0sq$(nr0sq)_h$(h).txt")
    open(out_path, "w") do io
        println(io, "# Move_all Δτ Sweep — bilayer")
        println(io, "# N=$(N), nr0sq=$(nr0sq), h=$(h), R0_opt=$(R0_opt), R_match=$(R_match)")
        println(io, "# delta_opt=$(delta_opt_moveall), Δτ_DMC=$(Δτ_DMC)")
        println(io, "#")
        println(io, "# Δτ\tAcceptance")
        for (δ, acc) in zip(delta_vals, acceptance_vals)
            println(io, "$(δ)\t$(acc)")
        end
        println(io, "#")
        println(io, "# Δτ_opt\tΔτ_DMC")
        println(io, "$(delta_opt_moveall)\t$(Δτ_DMC)")
    end
    println("✓ Saved to: $out_path")
end