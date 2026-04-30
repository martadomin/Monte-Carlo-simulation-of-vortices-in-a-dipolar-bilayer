# main.jl

using Random, Bessels, Plots, LaTeXStrings, ProgressMeter, Base.Threads, DelimitedFiles

# Source files
include(normpath(joinpath(@__DIR__, "src", "utils.jl")))
include(normpath(joinpath(@__DIR__, "src", "jastrow.jl")))
include(normpath(joinpath(@__DIR__, "src", "energy.jl")))
include(normpath(joinpath(@__DIR__, "src", "metropolis.jl")))
include(normpath(joinpath(@__DIR__, "src", "optimization.jl")))
include(normpath(joinpath(@__DIR__, "src", "observables.jl")))

# Parameters
num_part = 30
nr0_sq_vals = [16.0]  # main values for production runs
# nr0_sq_vals = [512.0, 768.0, 1024.0]  # for testing´

num_steps_coarse = 10^6
num_steps_fine = 10^7
num_steps_production = 10^7


println("Starting simulation...")

@threads for nr0_sq_val in nr0_sq_vals
    global nr0_sq = nr0_sq_val
    global L = sqrt(num_part / nr0_sq)
    global R_opt  # will be set by optimize_Rmatch.jl

    println("\n=== Simulation for N = $num_part, nr0^2 = $nr0_sq, L = $L ===")
    @time begin
        # First step: optimize R_match
        println("\n=== Stage 1: R_match Optimization ===")
        include(normpath(joinpath(@__DIR__, "scripts", "optimize_Rmatch.jl")))

        # Second step: production VMC run
        println("\n=== Stage 2: Production VMC Run ===")
        include(normpath(joinpath(@__DIR__, "scripts", "run_vmc.jl")))

        # Third step: plots
        println("\n=== Stage 3: Generating Plots ===")
        include(normpath(joinpath(@__DIR__, "scripts", "plot_results.jl")))
    end
end
println("N = $num_part, nr0^2 = $nr0_sq, L = $L")



println("\nTotal simulation completed.")