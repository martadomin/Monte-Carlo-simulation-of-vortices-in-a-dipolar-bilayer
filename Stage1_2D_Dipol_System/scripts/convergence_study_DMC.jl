# scripts/run_dmc.jl
using DelimitedFiles, Plots, LaTeXStrings, Statistics, ProgressMeter, Base.Threads

include(normpath(joinpath(@__DIR__, "..", "src", "utils.jl")))
include(normpath(joinpath(@__DIR__, "..", "src", "jastrow.jl")))
include(normpath(joinpath(@__DIR__, "..", "src", "energy.jl")))
include(normpath(joinpath(@__DIR__, "..", "src", "metropolis.jl")))
include(normpath(joinpath(@__DIR__, "..", "src", "dmc.jl")))
include(normpath(joinpath(@__DIR__, "..", "src", "observables.jl")))

# ── Parameters ────────────────────────────────────────────────────────────────
num_part      = 30
num_walkers   = 400
quadratic     = true
type_dmc      = quadratic ? "quadratic" : "linear"
total_time    = 10.0
equil_time    = 2.0     # fixed equilibration imaginary time
τ_factors = [1.0]  # factors to multiply Δτ_ref by

# ── Load reference Δτ table ───────────────────────────────────────────────────
ref_data     = readdlm(joinpath(@__DIR__, "..", "data", "sweep_results",
               "delta_optimal_N$(num_part)_DMC.txt"), '\t', Float64, skipstart=1)
println(joinpath(@__DIR__, "..", "data", "sweep_results",
               "delta_optimal_N$(num_part)_DMC.txt"))
nr0_sq_vals  = ref_data[:, 1]
print(nr0_sq_vals)
Δτ_ref_vals  = ref_data[:, 4]
print(Δτ_ref_vals)

# ── Main loop ─────────────────────────────────────────────────────────────────
for (nr0_sq, Δτ_ref) in zip(nr0_sq_vals, Δτ_ref_vals)
    
    L = sqrt(num_part / nr0_sq)

    # Load VMC data
    vmc_data      = readdlm(joinpath(@__DIR__, "..", "data", "results", "VMC",
                    "vmc_N$(num_part)_nr0sq$(nr0_sq).txt"), '\t', String, skipstart=1)
    R_opt         = parse(Float64, vmc_data[1, 4])
    E_ref_initial = parse(Float64, vmc_data[1, 5])
    Constants     = calculate_constants(L, R_opt)

    config_data = readdlm(joinpath(@__DIR__, "..", "data", "results", "VMC",
                  "vmc_config_N$(num_part)_nr0sq$(nr0_sq).txt"), '\t', Float64, skipstart=1)
    x_init = config_data[:, 1]
    y_init = config_data[:, 2]

    results_path = joinpath(@__DIR__, "..", "data", "results", "DMC",
                   "dmc_N$(num_part)_nr0sq$(nr0_sq)_Nw$(num_walkers)_$(type_dmc).txt")

    open(results_path, "w") do io
        println(io, "num_walkers\tΔτ\tE_dmc\tError_E_dmc")

        for factor in τ_factors
            Δτ         = Δτ_ref * factor
            num_steps  = max(10^4, round(Int, total_time / Δτ))
            num_equil  = max(2000, round(Int, equil_time / Δτ))

            println("\n═══ nr0_sq=$nr0_sq | Δτ=$(round(Δτ, sigdigits=3)) | steps=$num_steps | number of walkers=$num_walkers ═══")

            E_dmc, E_dmc_err, E_history = dmc(
                x_init, y_init,
                num_walkers, num_part, num_steps,
                Δτ, L, R_opt, Constants,
                E_ref_initial, num_walkers;
                num_equil   = num_equil,
                quadratic   = quadratic,
                plot_energy = false
            )

            # Block averaging
            block_sizes = [10, 20, 30, 40, 50, 100, 150, 200,
                           300, 400, 500, 600, 700, 800, 900, 1000,
                           1100, 1200, 1300, 1400, 1500, 1600, 1700,
                           1800, 1900, 2000]
            sigmas = [blocking_statistics(E_history, B)[2] for B in block_sizes]
            
            plateau  = detect_plateau(block_sizes, sigmas, window_size=4, rtol=0.02)

            p = plot(block_sizes, sigmas, marker=:o,
                     xlabel="Block Size", ylabel="Std. Error", title="Blocking Analysis (Δτ=$(round(Δτ, sigdigits=3)))")
            vline!([plateau], label="Plateau Block Size = $plateau", linestyle=:dash, color=:red)
            display(p)

            avg_E, σ = blocking_statistics(E_history, plateau)

            println("E/N = $(round(avg_E/num_part, digits=5)) ± $(round(σ/num_part, digits=5))")

            println(io, "$(num_walkers)\t$(Δτ)\t$(avg_E)\t$(σ)")
            flush(io)
        end
    end
    println("\nSaved: $results_path")
end

println("\nAll DMC runs completed.")