# Stage1_2D_Dipol_System/main_VMC.jl

include(joinpath(@__DIR__, "..", "common", "src", "CommonCore.jl"))
include(joinpath(@__DIR__, "src", "wavefunction.jl"))
include(joinpath(@__DIR__, "src", "optimization.jl"))
include(joinpath(@__DIR__, "config.jl"))

## ==== VARIATIONAL MONTE CARLO SIMULATION OF 2D DIPOLAR SYSTEM ==== ##

num_steps_coarse     = 10^6
num_steps_fine       = 10^6
num_steps_production = 10^7

plot_Rmatch_sweep_diagnostic = true   # set true to also show the R_match sweep plot per nr0_sq

stage_dir = joinpath(@__DIR__)

println("Starting simulation...")

for nr0_sq_val in nr0_sq_vals
    global nr0_sq = nr0_sq_val
    global L = sqrt(num_part / nr0_sq_val)

    println("\n=== Simulation for N = $num_part, nr0^2 = $nr0_sq_val, L = $L ===")

    production_path = result_path(stage_dir, "VMC", (N=num_part, nr0sq=nr0_sq))

    if skip_if_exists && isfile(production_path)
        println("Production VMC result already exists for N=$num_part, nr0²=$nr0_sq — skipping R_match optimization and production run entirely.")
        loaded_runs = load_run(production_path)
        avg_energy, sigma = combine_runs([r.result.avg_energy for r in loaded_runs], [r.result.sigma for r in loaded_runs])
        println("Combined E/N = ", avg_energy/num_part, " ± ", sigma/num_part, " (from $(length(loaded_runs)) run(s))")
    else
        @time begin
            println("\n=== Stage 1: R_match Optimization ===")
            include(joinpath(@__DIR__, "scripts", "optimize_Rmatch.jl"))

            plot_Rmatch_sweep_diagnostic && include(joinpath(@__DIR__, "scripts", "plot_Rmatch_sweep.jl"))

            println("\n=== Stage 2: Production VMC Run ===")
            include(joinpath(@__DIR__, "scripts", "run_vmc.jl"))
        end
    end
end

println("\n=== Final plot ===")
include(joinpath(@__DIR__, "scripts", "final_plot_VMC.jl"))

println("\nTotal simulation completed.")