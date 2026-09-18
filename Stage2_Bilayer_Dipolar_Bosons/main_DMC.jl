# Stage2_Bilayer_Dipolar_Bosons/main_DMC.jl

include(joinpath(@__DIR__, "..", "common", "src", "CommonCore.jl"))
include(joinpath(@__DIR__, "src", "wavefunction.jl"))
include(joinpath(@__DIR__, "src", "optimization.jl"))
include(joinpath(@__DIR__, "src", "h_sweep_analysis.jl"))
include(joinpath(@__DIR__, "config.jl"))

## ==== DMC SIMULATION OF 2D BILAYER DIPOLAR SYSTEM ==== ##
# Requires main_VMC.jl already run for every (nr0_sq, h) below — this
# pipeline loads R0_opt and the VMC reference energy from Stage2's own
# saved VMC results, and R_match from Stage1's, rather than recomputing.

quadratic = true
num_steps_dmc = 10^5

num_walkers_dtau_study = 200
num_walkers_vals = [20, 30, 40, 50, 75, 150, 200, 300, 400]

println("Starting DMC simulations...")

for nr0_sq in nr0_sq_vals
    global L = sqrt(num_part / nr0_sq)
    println("\n=== nr0^2 = $nr0_sq, L = $L ===")

    for h in h_vals
        global h
        println("\n--- h = $h ---")
        @time begin
            println("\n--- Step 1: Tune delta for DMC ---")
            include(joinpath(@__DIR__, "scripts", "tune_delta_dmc.jl"))

            println("\n--- Step 2: Δτ convergence ---")
            include(joinpath(@__DIR__, "scripts", "dtau_convergence.jl"))

            println("\n--- Step 3: Walker-count convergence (final result) ---")
            include(joinpath(@__DIR__, "scripts", "num_walkers_convergence.jl"))
        end
    end

    println("\n=== h-dependence DMC plots for nr0^2 = $nr0_sq ===")
    include(joinpath(@__DIR__, "scripts", "final_inset_DMC.jl"))
    include(joinpath(@__DIR__, "scripts", "final_fig1_DMC.jl"))
end

println("\nTotal DMC simulation completed.")