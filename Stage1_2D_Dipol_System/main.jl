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
num_part             = 30
num_steps_coarse     = 10^5
num_steps_fine       = 10^6
num_steps_production = 10^6

nr0_sq_values = [16.0, 32.0, 48.0, 64.0, 96.0, 128.0, 196.0, 256.0, 384.0, 512.0, 768.0, 1024.0]

global nr0_sq  # declare before the loop
global L
global R_opt   # any other variables used in the scripts

println("Starting simulation...")
println("N = $num_part")

for nr0_sq_val in nr0_sq_values
    global nr0_sq = nr0_sq_val
    global L = sqrt(num_part / nr0_sq)

    println("\n========================================")
    println("nr0^2 = $nr0_sq, L = $L")
    println("========================================")

    @time begin
        println("\n=== Stage 1: R_match Optimization ===")
        include(normpath(joinpath(@__DIR__, "scripts", "optimize_Rmatch.jl")))
        println("\n=== Stage 2: Production VMC Run ===")
        include(normpath(joinpath(@__DIR__, "scripts", "run_vmc.jl")))
        println("\n=== Stage 3: Generating Plots ===")
        include(normpath(joinpath(@__DIR__, "scripts", "plot_results.jl")))
    end
    println("Completed nr0^2 = $nr0_sq")
end

println("\nAll simulations completed.")