# Stage1_2D_Dipol_System/main_DMC.jl

include(joinpath(@__DIR__, "..", "common", "src", "CommonCore.jl"))
include(joinpath(@__DIR__, "src", "wavefunction.jl"))
include(joinpath(@__DIR__, "src", "optimization.jl"))
include(joinpath(@__DIR__, "config.jl"))
## ==== DIFFUSION MONTE CARLO SIMULATION OF 2D DIPOLAR SYSTEM ==== ##

quadratic = true
num_steps_dmc = 10^5

num_walkers_dtau_study = 200                       # fixed, for the Δτ study only
num_walkers_vals = [20, 30, 40, 50, 75, 150, 200, 300, 400]  # swept for the final extrapolation

stage_dir = joinpath(@__DIR__)

println("Starting DMC simulations...")

for nr0_sq in nr0_sq_vals
    global L = sqrt(num_part / nr0_sq)

    # R_opt for this nr0_sq must already exist — from having run main_VMC.jl
    # first, which saves it via optimize_Rmatch.jl's own save_run call.
    R_opt_run = load_run(result_path(stage_dir, "Rmatch_optimum", (N=num_part, L=L)); run=1)
    global R_opt = R_opt_run.result.R_opt

    println("\n=== DMC Simulation for N = $num_part, nr0^2 = $nr0_sq, L = $L, R_opt = $R_opt ===")
    @time begin
        println("\n--- Step 1: Tune delta for DMC ---")
        include(joinpath(@__DIR__, "scripts", "tune_delta_dmc.jl"))

        println("\n--- Step 2: Δτ convergence ---")
        include(joinpath(@__DIR__, "scripts", "dtau_convergence.jl"))

        println("\n--- Step 3: Walker-count convergence (final result) ---")
        include(joinpath(@__DIR__, "scripts", "num_walkers_convergence.jl"))
    end
end

println("\n=== Final plot ===")
include(joinpath(@__DIR__, "scripts", "final_plot_DMC.jl"))

println("\nTotal DMC simulation completed.")