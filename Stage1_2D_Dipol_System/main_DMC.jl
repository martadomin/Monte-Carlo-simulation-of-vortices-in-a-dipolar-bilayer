using DelimitedFiles, Plots, LaTeXStrings, Statistics, ProgressMeter, Base.Threads

# Source files
include(normpath(joinpath(@__DIR__, "src", "utils.jl")))
include(normpath(joinpath(@__DIR__, "src", "jastrow.jl")))
include(normpath(joinpath(@__DIR__, "src", "energy.jl")))
include(normpath(joinpath(@__DIR__, "src", "metropolis.jl")))
include(normpath(joinpath(@__DIR__, "src", "optimization.jl")))
include(normpath(joinpath(@__DIR__, "src", "observables.jl")))
include(normpath(joinpath(@__DIR__, "src", "dmc.jl")))

# Parameters
num_part = 30
# nr0_sq_vals = [256.0, 384.0, 512.0, 768.0, 1024.0]
nr0_sq_vals = [16.0, 32.0, 48.0, 64.0, 96.0, 128.0]

num_steps_MC = 10^5
num_walkers_vals = [1000]
Δτ_vals = [1e-4]

println("Starting DMC simulations...")

for nr0_sq_val in nr0_sq_vals
    global nr0_sq = nr0_sq_val
    global L = sqrt(num_part / nr0_sq)
    global R_opt
    global nr0_sq_32 = num_part * nr0_sq^(3/2)

    println("\n=== DMC Simulation for N = $num_part, nr0^2 = $nr0_sq, L = $L ===")
    include(normpath(joinpath(@__DIR__, "scripts", "run_dmc.jl")))
end

include(normpath(joinpath(@__DIR__, "scripts", "plot_results_DMC.jl")))

