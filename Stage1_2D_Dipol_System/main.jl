# main.jl

# Dependencies
using Random, Bessels, Plots, LaTeXStrings, ProgressMeter, Base.Threads, DelimitedFiles

# Source files
include(normpath(joinpath(@__DIR__, "src", "utils.jl")))
include(normpath(joinpath(@__DIR__, "src", "jastrow.jl")))
include(normpath(joinpath(@__DIR__, "src", "energy.jl")))
include(normpath(joinpath(@__DIR__, "src", "metropolis.jl")))
include(normpath(joinpath(@__DIR__, "src", "optimization.jl")))

# Parameters
num_part             = 30
nr0_sq               = 16.0
L                    = sqrt(num_part / nr0_sq)
num_steps_coarse     = 10^5
num_steps_fine       = 10^6
num_steps_production = 10^7

# First step: optimize R_match
include(normpath(joinpath(@__DIR__, "scripts", "optimize_Rmatch.jl")))

# Second step: production VMC run
include(normpath(joinpath(@__DIR__, "scripts", "run_vmc.jl")))

# Third step: plots
include(normpath(joinpath(@__DIR__, "scripts", "plot_results.jl")))