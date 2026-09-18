# Stage2_Bilayer_Dipolar_Bosons/main_VMC.jl

include(joinpath(@__DIR__, "..", "common", "src", "CommonCore.jl"))
include(joinpath(@__DIR__, "src", "wavefunction.jl"))
include(joinpath(@__DIR__, "src", "optimization.jl"))
include(joinpath(@__DIR__, "src", "h_sweep_analysis.jl"))
include(joinpath(@__DIR__, "config.jl"))

## ==== VMC SIMULATION OF 2D BILAYER DIPOLAR SYSTEM ==== ##
num_steps_coarse     = 10^5
num_steps_fine       = 10^6
num_steps_production = 10^6

plot_R0_sweep_diagnostic = false

println("Starting simulation...")

for nr0_sq in nr0_sq_vals
    global L = sqrt(num_part / nr0_sq)
    println("\n=== nr0^2 = $nr0_sq, L = $L ===")

    for h in h_vals
        global h
        println("\n--- h = $h ---")
        @time begin
            println("\n=== Stage 1: R0 Optimization ===")
            include(joinpath(@__DIR__, "scripts", "optimize_R0.jl"))

            plot_R0_sweep_diagnostic && include(joinpath(@__DIR__, "scripts", "plot_R0_sweep.jl"))

            println("\n=== Stage 2: Production VMC Run ===")
            include(joinpath(@__DIR__, "scripts", "run_vmc.jl"))
        end
    end

    println("\n=== h-dependence plots for nr0^2 = $nr0_sq ===")
    include(joinpath(@__DIR__, "scripts", "final_inset_VMC.jl"))
    include(joinpath(@__DIR__, "scripts", "final_fig1_VMC.jl"))
end

println("\nTotal simulation completed.")