# Stage3_Vortex_Excitations/main_VMC.jl

include(joinpath(@__DIR__, "..", "common", "src", "CommonCore.jl"))
include(joinpath(@__DIR__, "src", "wavefunction.jl"))
include(joinpath(@__DIR__, "src", "loaders.jl"))
include(joinpath(@__DIR__, "config.jl"))

## ==== VMC SIMULATION OF THE VORTEX BILAYER — OFFSET SWEEP ==== ##
# Requires main_VMC.jl already run for Stage1 (this N_half, L) and
# Stage2 (this num_part, nr0_sq, h) — no optimization here.

println("Starting vortex offset sweep (VMC)...")

for d in d_vals
    global d
    @time include(joinpath(@__DIR__, "scripts", "run_vmc.jl"))
end

println("\n=== Final plot: E_vortex(d)/N ===")
include(joinpath(@__DIR__, "scripts", "final_plot_vortex_VMC.jl"))

println("\nTotal simulation completed.")