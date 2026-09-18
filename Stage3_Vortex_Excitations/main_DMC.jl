# Stage3_Vortex_Excitations/main_DMC.jl

include(joinpath(@__DIR__, "..", "common", "src", "CommonCore.jl"))
include(joinpath(@__DIR__, "src", "wavefunction.jl"))
include(joinpath(@__DIR__, "src", "loaders.jl"))
include(joinpath(@__DIR__, "config.jl"))

## ==== DMC SIMULATION OF THE VORTEX BILAYER — OFFSET SWEEP ==== ##
# Requires main_VMC.jl already run for this Stage3 sweep (VMC_vortex,
# per d) — E_ref_initial is loaded from there.

println("Starting vortex offset sweep (DMC)...")

for d in d_vals
    global d
    @time begin
        println("\n--- Step 1: Tune delta for DMC, d=$d ---")
        include(joinpath(@__DIR__, "scripts", "tune_delta_dmc.jl"))

        println("\n--- Step 2: Δτ convergence, d=$d ---")
        include(joinpath(@__DIR__, "scripts", "dtau_convergence.jl"))

        println("\n--- Step 3: Walker-count convergence (final result), d=$d ---")
        include(joinpath(@__DIR__, "scripts", "num_walkers_convergence.jl"))
    end
end

println("\n=== Final plot: E_vortex(d)/N ===")
include(joinpath(@__DIR__, "scripts", "final_plot_vortex_DMC.jl"))

println("\nTotal DMC simulation completed.")