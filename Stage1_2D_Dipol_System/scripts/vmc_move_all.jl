# scripts/run_vmc.jl
include(normpath(joinpath(@__DIR__, "..", "src", "utils.jl")))
include(normpath(joinpath(@__DIR__, "..", "src", "jastrow.jl")))
include(normpath(joinpath(@__DIR__, "..", "src", "energy.jl")))
include(normpath(joinpath(@__DIR__, "..", "src", "metropolis.jl")))
include(normpath(joinpath(@__DIR__, "..", "src", "observables.jl")))

using Plots, LaTeXStrings, Base.Threads, DelimitedFiles, Statistics


nr0_sq_vals = [96.0]
num_part = 30

results_path = joinpath(@__DIR__, "..", "data", "sweep_results",
               "delta_optimal_N$(num_part)_DMC.txt")
open(results_path, "w") do io
    println(io, "nr0_sq\tdelta_opt\tacceptance\tDelta_tau_DMC")
end

for nr0_sq in nr0_sq_vals
    L= sqrt(num_part / nr0_sq)
    nr0_sq_32 = num_part * nr0_sq^(3/2)

    # ── Load R_opt ────────────────────────────────────────────────────────────────
    rmatch_path = joinpath(@__DIR__, "..", "data", "sweep_results",
                "Rmatch_sweep_N$(num_part)_nr0sq$(nr0_sq).txt")

    R_opt_ref = Ref(0.0)
    open(rmatch_path, "r") do io
        section = ""
        for line in eachline(io)
            startswith(line, "#") && (section = line; continue)
            isempty(strip(line)) && continue
            line == "R_opt\tE_opt" && continue
            if occursin("Optimal R_match", section)
                vals = parse.(Float64, split(line, "\t"))
                R_opt_ref[] = vals[1]
            end
        end
    end
    R_opt = R_opt_ref[]
    println("Loaded R_opt = $(R_opt)")

    Constants = calculate_constants(L, R_opt)
    x_coord, y_coord = random_initial_config(num_part, L, "Uniform")

    # ── Delta sweep ───────────────────────────────────────────────────────────────
    println("\n--- Sweeping delta for 50% acceptance (move_all) ---")

    delta_vals = [0.003, 0.0031, 0.0032, 0.0033, 0.0034, 0.0035, 0.0036, 0.0037, 0.0038, 0.0039,
                  0.004]

    acceptance_vals = Float64[]
    num_sweep_steps = 10^4

    for δ in delta_vals
        _, _, _, _, _, _, _, _, _, acc, _, _, _, _ = metropolis(
            num_part, num_sweep_steps, δ, L, R_opt, Constants;
            x_init       = x_coord,
            y_init       = y_coord,
            final_energy_plot = false,
            progress     = true,
            move_all     = true
        )
        push!(acceptance_vals, acc)
        println("δ = $(rpad(δ, 6)) → Acceptance = $(round(acc*100, digits=3))%")
    end

    # ── Find delta closest to 50% ─────────────────────────────────────────────────
    best_idx = argmin(abs.(acceptance_vals .- 0.5))
    delta_opt = delta_vals[best_idx]
    Δτ_DMC = delta_opt^2 / (2 * 0.5)  # D = ħ²/(2m) with ħ=1, m=1
    println("\nOptimal delta = $(delta_opt)  (acceptance = $(round(acceptance_vals[best_idx]*100, digits=3))%)")
    println("\nFor the DMC run: δ = √{DΔτ} → Δτ = δ^2/D = " * string(round(Δτ_DMC*0.1, digits=7)))

    # ── Plot acceptance vs delta ──────────────────────────────────────────────────
    p_sweep = plot(delta_vals, acceptance_vals;
        xlabel    = L"\delta",
        ylabel    = "Acceptance ratio",
        title     = "Delta sweep (move all), N=$(num_part), nr₀²=$(nr0_sq)",
        marker    = :circle,
        label     = "acceptance",
        xscale    = :log10,
        legend    = :topright
    )
    hline!(p_sweep, [0.5]; linestyle=:dash, color=:red, label="50%")
    vline!(p_sweep, [delta_opt]; linestyle=:dot, color=:black, label="optimal δ")
    display(p_sweep)

    # ── Save sweep results ─────────────────────────────────────────────────────────
    open(results_path, "a") do io
        println(io, "$(nr0_sq)\t$(delta_opt)\t$(acceptance_vals[best_idx])\t$(Δτ_DMC)")
    end
end

# # ── Production VMC ────────────────────────────────────────────────────────────
# println("\n--- Running production VMC ($num_steps_production steps, delta=$(delta_opt)) ---")

# energies_vmc, energies_drift_vmc, energies_laplacian_vmc,
# E_tot, _, E_drift, E_laplacian, _, _, acceptance_ratio,
# r_vals, g_r_normalized, x_final, y_final = metropolis(
#     num_part, num_steps_production, delta_opt, L, R_opt, Constants;
#     x_init            = x_coord,
#     y_init            = y_coord,
#     final_energy_plot = false,
#     plot_every        = 10^3,
#     progress          = true,
#     num_bins          = 100,
#     move_all          = true
# )

# # ── Block averaging ───────────────────────────────────────────────────────────
# block_sizes = [10, 20, 30, 40, 50, 100, 150, 200, 300, 400, 500,
#                600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400,
#                1500, 1600, 1700, 1800, 1900, 2000]
# sigmas            = Float64[]
# sigmas_drift      = Float64[]
# sigmas_laplacian  = Float64[]

# for B in block_sizes
#     _, σ  = blocking_statistics(energies_vmc, B)
#     push!(sigmas, σ)
# end

# plateau_std = detect_plateau(block_sizes, sigmas, window_size=4, rtol=0.05)

# avg_energy, sigma = blocking_statistics(energies_vmc, plateau_std)

# println("\n--- Results ---")
# println("E/N/(nr0^2)^(3/2)          = $(avg_energy/nr0_sq_32) ± $(sigma/nr0_sq_32)")
# println("Acceptance ratio           = $(acceptance_ratio)")
# println("Optimal delta              = $(delta_opt)")

# # ── Save results ──────────────────────────────────────────────────────────────
# results_path = joinpath(@__DIR__, "..", "data", "results", "VMC",
#                "vmc_N$(num_part)_nr0sq$(nr0_sq).txt")
# open(results_path, "w") do io
#     println(io, "num_part\tnr0_sq\tL\tR_opt\tE_tot\tError")
#     println(io, "$(num_part)\t$(nr0_sq)\t$(L)\t$(Float64(R_opt))\t$(avg_energy)\t$(sigma)")
# end

# config_path = joinpath(@__DIR__, "..", "data", "results", "VMC",
#               "vmc_config_N$(num_part)_nr0sq$(nr0_sq).txt")
# open(config_path, "w") do io
#     println(io, "x\ty")
#     for (x, y) in zip(x_final, y_final)
#         println(io, "$(x)\t$(y)")
#     end
# end

# gr_path = joinpath(@__DIR__, "..", "data", "results", "VMC",
#           "gr_N$(num_part)_nr0sq$(nr0_sq).txt")
# open(gr_path, "w") do io
#     println(io, "r\tg(r)")
#     for (r, g) in zip(r_vals, g_r_normalized)
#         println(io, "$(r)\t$(g)")
#     end
# end
# println("Results saved.")