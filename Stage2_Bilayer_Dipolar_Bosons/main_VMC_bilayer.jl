"""
main_VMC_Stage2.jl

Entry point for Stage II bilayer VMC simulations.
Sets all physical and numerical parameters, then for each (N, h) pair:
  1. Runs the R0 sweep (optimize_R0.jl) to find R0_opt
  2. Runs the production VMC (run_vmc_bilayer.jl) at R0_opt

No parameters are hardcoded inside the subscripts — everything flows
from the globals defined here.
"""

using Random, Plots, LaTeXStrings, Printf, Base.Threads

include(normpath(joinpath(@__DIR__, "src", "utils.jl")))
include(normpath(joinpath(@__DIR__, "src", "jastrow.jl")))
include(normpath(joinpath(@__DIR__, "src", "energy.jl")))
include(normpath(joinpath(@__DIR__, "src", "metropolis.jl")))
include(normpath(joinpath(@__DIR__, "src", "shooting_method.jl")))
include(normpath(joinpath(@__DIR__, "src", "observables.jl")))

"""
    read_optimal_stage1(path)

Parse the trailing "# Optimal R_match" block of a Stage I Rmatch sweep file:
    # Optimal R_match
    R_opt	E_opt
    <value>	<value>
Returns (R_opt, E_opt_total) as Float64. E_opt_total is the TOTAL energy
of the N÷2 single-layer system (not per-particle).
"""
function read_optimal_stage1(path::String)
    lines = readlines(path)
    idx = findfirst(l -> occursin("Optimal R_match", l), lines)
    idx === nothing && error("'# Optimal R_match' block not found in $path")
    vals = split(strip(lines[idx + 2]), '\t')
    return parse(Float64, vals[1]), parse(Float64, vals[2])
end


println("\n" * "="^70)
println("STAGE 2: BILAYER DIPOLAR VMC")
println("="^70 * "\n")

# ──────────────────────────────────────────────────────────────────
# PHYSICAL PARAMETERS
# ──────────────────────────────────────────────────────────────────

global N       = 60               # Total particles (N/2 per layer)
global nr0sq   = 1.0              # Areal density × r₀²
global L = sqrt(N / (nr0sq))

# Load R_match from Stage I sweep results (single-layer, N/2 particles)
rmatch_path = joinpath(@__DIR__, "..", "Stage1_2D_Dipol_System", "data", "sweep_results",
                       "Rmatch_sweep_N$(N÷2)_nr0sq$(nr0sq/2).txt")
                       

global R_match, E_opt_single_layer_total = read_optimal_stage1(rmatch_path)
println("✓ R_match = $R_match r₀  (loaded from $(basename(rmatch_path)))")

# ──────────────────────────────────────────────────────────────────
# SHOOTING METHOD PARAMETERS (fixed)
# ──────────────────────────────────────────────────────────────────

global r_min     = 1e-6
global Δ_shoot   = 1e-4
global tol_shoot = 1e-10

# ──────────────────────────────────────────────────────────────────
# VMC PARAMETERS
# ──────────────────────────────────────────────────────────────────

global num_steps_coarse    = 10^6  # MC steps per R0 — coarse sweep
global num_steps_fine      = 10^6 # MC steps per R0 — fine sweep
global num_steps_production = 10^6  # MC steps — production run
global num_tune_steps      = 10000  # Steps for delta tuning
global n_points_sweep      = 20    # Points per sweep stage

# Interlayer separations to sweep
global h_vals = range(0.3, 1.5, step = 0.1)

println(collect(h_vals))
# ──────────────────────────────────────────────────────────────────
# EXACT DIMER BINDING ENERGIES
# ──────────────────────────────────────────────────────────────────

include(normpath(joinpath(@__DIR__, "scripts", "find_binding_energy.jl")))

global Constants = calculate_constants(L, R_match)

# ──────────────────────────────────────────────────────────────────
# DERIVED QUANTITIES (recomputed inside loop for each h)
# ──────────────────────────────────────────────────────────────────

# Storage for results across h values
global results = Dict{Float64, NamedTuple}()

for h_val in h_vals
    global h = h_val

    println("\n" * "="^70)
    println("h = $h r₀  |  N = $N  |  nr0sq = $nr0sq  |  L = $(round(L, digits=4)) r₀")
    println("="^70)

    # ── Step 1: R0 optimization ──────────────────────────────────
    println("\n--- R0 Optimization ---")
    include(normpath(joinpath(@__DIR__, "scripts", "optimize_R0.jl")))
    # optimize_R0.jl sets: R0_opt, E_opt, err_opt, energy_b_opt

    # Plot the R0 sweep results for this h
    include(normpath(joinpath(@__DIR__, "scripts", "plot_R0_sweep.jl")))

    println("\n  → R0_opt = $(round(R0_opt, digits=6)) r₀")
    println("  → E/N    = $(round(E_opt,  digits=6)) ± $(round(err_opt, digits=6)) ħ²/mr₀²")

    # ── Step 2: Production VMC at R0_opt ─────────────────────────
    println("\n--- Production VMC ---")
    prod_results = include(normpath(joinpath(@__DIR__, "scripts", "run_VMC_bilayer.jl")))
    E_production   = prod_results.E_production
    err_production  = prod_results.err_production

    # ── Store results ─────────────────────────────────────────────
    results[h_val] = (
        R0_opt        = R0_opt,
        energy_b_opt  = energy_b_opt,
        E_vmc         = E_production,
        err_vmc       = err_production,
    )

    println("\n✓ h = $h done:  E/N = $(round(E_production, digits=6)) ± $(round(err_production, digits=6))")
end

# ──────────────────────────────────────────────────────────────────
# SUMMARY
# ──────────────────────────────────────────────────────────────────

println("\n" * "="^70)
println("SUMMARY — nr0sq = $nr0sq, N = $N")
println("="^70)
println("  h (r₀)   R0_opt (r₀)   ε_b          E/N ± err  (ħ²/mr₀²)")
println("-"^70)
for h_val in h_vals
    haskey(results, h_val) || continue
    r = results[h_val]
    @printf("  %-8.3f  %-12.6f  %-12.6f  %.6f ± %.6f\n",
            h_val, r.R0_opt, r.energy_b_opt, r.E_vmc, r.err_vmc)
end

# ──────────────────────────────────────────────────────────────────
# SAVE SUMMARY
# ──────────────────────────────────────────────────────────────────

summary_path = joinpath(@__DIR__, "data", "VMC_Stage2_N$(N)_nr0sq$(nr0sq)_summary.txt")
mkpath(dirname(summary_path))

open(summary_path, "w") do io
    println(io, "# Stage 2 VMC summary")
    println(io, "# N=$(N), nr0sq=$(nr0sq), R_match=$(R_match)")
    println(io, "# h\tR0_opt\teps_b\tE_vmc\terr_vmc")
    for h_val in h_vals
        haskey(results, h_val) || continue
        r = results[h_val]
        println(io, "$(h_val)\t$(r.R0_opt)\t$(r.energy_b_opt)\t$(r.E_vmc)\t$(r.err_vmc)")
    end
end

println("\n✓ Summary saved to: $summary_path")
println("="^70 * "\n")

# ──────────────────────────────────────────────────────────────────
# PLOT RESULTS
# ──────────────────────────────────────────────────────────────────

# Figure inset
include(normpath(joinpath(@__DIR__, "scripts", "plot_inset_fig1.jl")))

# Figure 1
include(normpath(joinpath(@__DIR__, "scripts", "plot_fig1.jl")))