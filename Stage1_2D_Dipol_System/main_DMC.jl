# Stage1_2D_Dipol_System/main_DMC.jl

include(joinpath(@__DIR__, "..", "common", "src", "CommonCore.jl"))
include(joinpath(@__DIR__, "src", "wavefunction.jl"))
include(joinpath(@__DIR__, "src", "optimization.jl"))
include(joinpath(@__DIR__, "config.jl"))

## ==== DMC SIMULATION OF 2D DIPOLAR SYSTEM ==== ##
# Requires main_VMC.jl already run for every nr0_sq below — this pipeline
# loads R_opt and the VMC reference energy from Stage1's own saved VMC
# results, rather than recomputing.

quadratic = false
num_steps_dmc = 10^5

num_walkers_dtau_study = 200
num_walkers_vals = [20, 30, 40, 50, 75, 150, 200, 300, 400]

stage_dir = joinpath(@__DIR__)

println("Starting DMC simulations...")

for nr0_sq_val in nr0_sq_vals
    global nr0_sq = nr0_sq_val
    global L = sqrt(num_part / nr0_sq_val)
    println("\n=== nr0^2 = $nr0_sq, L = $L ===")

    # Load R_opt from Stage1's own saved VMC optimization — must already
    # exist, from running main_VMC.jl for this num_part/L first.
    rmatch_path = result_path(stage_dir, "Rmatch_optimum", (N=num_part, L=L))
    if !isfile(rmatch_path)
        error("No saved R_match for N=$num_part, L=$L. Run main_VMC.jl first.")
    end
    global R_opt = load_run(rmatch_path; run=1).result.R_opt

    dmc_path = result_path(stage_dir, "DMC", (N=num_part, nr0sq=nr0_sq))

    if skip_if_exists && isfile(dmc_path)
        println("Final DMC result already exists for N=$num_part, nr0²=$nr0_sq — skipping.")
        loaded_runs = load_run(dmc_path)
        avg_E, sig = combine_runs([r.result.E_extrapolated for r in loaded_runs], [r.result.E_extrapolated_err for r in loaded_runs])
        println("Combined E_extrapolated = ", avg_E/num_part, " ± ", sig/num_part, " (from $(length(loaded_runs)) run(s))")
    else
        @time begin
            println("\n--- Step 1: Tune delta for DMC ---")
            include(joinpath(@__DIR__, "scripts", "tune_delta_dmc.jl"))

            println("\n--- Step 2: Δτ convergence ---")
            include(joinpath(@__DIR__, "scripts", "dtau_convergence.jl"))

            println("\n--- Step 3: Walker-count convergence (final result) ---")
            include(joinpath(@__DIR__, "scripts", "num_walkers_convergence.jl"))
        end
    end
end

println("\n=== Final plot ===")
include(joinpath(@__DIR__, "scripts", "final_plot_DMC.jl"))

println("\nTotal DMC simulation completed.")