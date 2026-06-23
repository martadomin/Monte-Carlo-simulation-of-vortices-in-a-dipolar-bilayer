using DelimitedFiles, Plots, LaTeXStrings, Statistics, ProgressMeter, Base.Threads

# Source files
include(normpath(joinpath(@__DIR__, "src", "utils.jl")))
include(normpath(joinpath(@__DIR__, "src", "jastrow.jl")))
include(normpath(joinpath(@__DIR__, "src", "energy.jl")))
include(normpath(joinpath(@__DIR__, "src", "metropolis.jl")))
include(normpath(joinpath(@__DIR__, "src", "optimization.jl")))
include(normpath(joinpath(@__DIR__, "src", "observables.jl")))
include(normpath(joinpath(@__DIR__, "src", "dmc.jl")))

# ── Parameters ────────────────────────────────────────────────────────────────
num_part      = 30
num_walkers   = 500
quadratic     = true
type_dmc      = quadratic ? "quadratic" : "linear"
total_time    = 10.0
equil_time    = 2.0     # fixed equilibration imaginary time
τ_factors = [1.0]  # factors to multiply Δτ_ref by

# ── Load reference Δτ table ───────────────────────────────────────────────────
file_path = joinpath(@__DIR__, "data", "sweep_results", "delta_optimal_N$(num_part)_DMC.txt")
if isfile(file_path)
    ref_data     = readdlm(file_path, '\t', Float64, skipstart=1)
    println(file_path)
    nr0_sq_vals  = ref_data[:, 1]
    print(nr0_sq_vals)
    Δτ_ref_vals  = ref_data[:, 4]
else
    Δτ = 10^-4
end

println("Starting DMC simulations...")

for nr0_sq_val in nr0_sq_vals
    global nr0_sq = nr0_sq_val
    global L = sqrt(num_part / nr0_sq)
    global R_opt
    global nr0_sq_32 = num_part * nr0_sq^(3/2)
    Δτ_base = Δτ_ref_vals[findfirst(x -> isapprox(x, nr0_sq; atol=1e-8), nr0_sq_vals)]
    global Δτ_vals = [Δτ_base * f for f in τ_factors]
    global num_walkers_vals = [num_walkers]
    global num_steps_dmc = max(10^4, round(Int, total_time / minimum(Δτ_vals)))

    println("\n=== DMC Simulation for N = $num_part, nr0^2 = $nr0_sq, L = $L ===")
    include(normpath(joinpath(@__DIR__, "scripts", "run_dmc.jl")))
end

# include(normpath(joinpath(@__DIR__, "scripts", "plot_results_DMC.jl")))

