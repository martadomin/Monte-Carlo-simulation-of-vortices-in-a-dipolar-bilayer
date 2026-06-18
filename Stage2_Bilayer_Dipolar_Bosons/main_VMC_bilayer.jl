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


println("\n" * "="^70)
println("STAGE 2: BILAYER DIPOLAR VMC")
println("="^70 * "\n")

# ──────────────────────────────────────────────────────────────────
# PHYSICAL PARAMETERS
# ──────────────────────────────────────────────────────────────────

global N       = 60               # Total particles (N/2 per layer)
global nr0sq   = 1.0              # Areal density × r₀²
global L = sqrt(N / nr0sq)

# Load R_match from Stage I sweep results (single-layer, N/2 particles)
rmatch_path = joinpath(@__DIR__, "..", "Stage1_2D_Dipol_System", "data", "sweep_results",
                       "Rmatch_sweep_N$(N÷2)_nr0sq$(nr0sq/2).txt")
print(rmatch_path)

global R_match = let
    last_line = ""
    for line in eachline(rmatch_path)
        if !startswith(line, "#") && !isempty(strip(line))
            last_line = line
        end
    end
    parse(Float64, split(last_line, "\t")[1])
end
println("✓ R_match = $R_match r₀  (loaded from $(basename(rmatch_path)))")

# Interlayer separations to sweep
global h_vals = range(0.3, 1.5, length=21)  # add h values here as the sweeps finish

# ──────────────────────────────────────────────────────────────────
# SHOOTING METHOD PARAMETERS (fixed)
# ──────────────────────────────────────────────────────────────────

global r_min     = 1e-6
global Δ_shoot   = 1e-4
global tol_shoot = 1e-10

# ──────────────────────────────────────────────────────────────────
# VMC PARAMETERS
# ──────────────────────────────────────────────────────────────────

global num_steps_coarse    = 10^5    # MC steps per R0 — coarse sweep
global num_steps_fine      = 3*10^5  # MC steps per R0 — fine sweep
global num_steps_production = 10^6   # MC steps — production run
global num_tune_steps      = 5000    # Steps for delta tuning
global n_points_sweep      = 12      # Points per sweep stage

# ──────────────────────────────────────────────────────────────────
# EXACT DIMER BINDING ENERGIES
# ──────────────────────────────────────────────────────────────────

include(normpath(joinpath(@__DIR__, "scripts", "find_binding_energy.jl")))

# ──────────────────────────────────────────────────────────────────
# DERIVED QUANTITIES (recomputed inside loop for each h)
# ──────────────────────────────────────────────────────────────────


# Storage for results across h values
results = Dict{Float64, NamedTuple}()

for h_val in h_vals
    global h = h_val
    global Constants = calculate_constants(L, R_match)

    println("\n" * "="^70)
    println("h = $h r₀  |  N = $N  |  nr0sq = $nr0sq  |  L = $(round(L, digits=4)) r₀")
    println("="^70)

    # ── Step 1: R0 optimization ──────────────────────────────────
    println("\n--- R0 Optimization ---")
    include(normpath(joinpath(@__DIR__, "scripts", "optimize_R0.jl")))
    # optimize_R0.jl sets: R0_opt, E_opt, err_opt, energy_b_opt

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