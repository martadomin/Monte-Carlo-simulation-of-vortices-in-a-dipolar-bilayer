# Stage1_2D_Dipol_System/main_VMC.jl

include(joinpath(@__DIR__, "..", "common", "src", "CommonCore.jl"))
include(joinpath(@__DIR__, "src", "wavefunction.jl"))
include(joinpath(@__DIR__, "src", "optimization.jl"))
include(joinpath(@__DIR__, "config.jl"))

## ==== VARIATIONAL MONTE CARLO SIMULATION OF 2D DIPOLAR SYSTEM ==== ##

num_steps_coarse     = 10^5
num_steps_fine       = 10^6
num_steps_production = 10^6

plot_Rmatch_sweep_diagnostic = false   # set true to also show the R_match sweep plot per nr0_sq

stage_dir = joinpath(@__DIR__)

println("Starting simulation...")

for nr0_sq in nr0_sq_vals
    global L = sqrt(num_part / nr0_sq)
    println("\n=== Simulation for N = $num_part, nr0^2 = $nr0_sq, L = $L ===")
    @time begin
        println("\n=== Stage 1: R_match Optimization ===")
        include(joinpath(@__DIR__, "scripts", "optimize_Rmatch.jl"))

        plot_Rmatch_sweep_diagnostic && include(joinpath(@__DIR__, "scripts", "plot_Rmatch_sweep.jl"))

        println("\n=== Stage 2: Production VMC Run ===")
        include(joinpath(@__DIR__, "scripts", "run_vmc.jl"))
    end
end

println("\n=== Final plot ===")
include(joinpath(@__DIR__, "scripts", "final_plot_VMC.jl"))

println("\nTotal simulation completed.")