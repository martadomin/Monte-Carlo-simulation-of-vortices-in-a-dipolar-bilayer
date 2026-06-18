"""
run_vmc_bilayer.jl

Production VMC run for the bilayer system at the optimal R0.
Reads globals from main_VMC_Stage2.jl:
  N, nr0sq, L, h, R_match, Constants, R0_opt, energy_b_opt,
  r_min, Δ_shoot, tol_shoot,
  num_steps_production, num_tune_steps

Sets on exit:
  E_production, err_production
"""

println("="^70)
println("PRODUCTION VMC  (h = $h r₀,  R0 = $(round(R0_opt, digits=6)) r₀)")
println("="^70 * "\n")

# Build fAB at R0_opt
r_grid_prod, psi_prod, u_prime_prod, u_doubleprime_prod, _ =
    build_fAB(h, R0_opt, r_min, Δ_shoot, tol_shoot, nr0sq, N)

println("fAB(r_min) = $(round(psi_prod[1],   digits=6))")
println("fAB(R0)    = $(round(psi_prod[end], digits=6))")

# Initial configuration
x_coord, y_coord = random_initial_config(N, L, "Uniform")
x_A = x_coord[1:N÷2];      y_A = y_coord[1:N÷2]
x_B = x_coord[N÷2+1:end];  y_B = y_coord[N÷2+1:end]

# Tune delta
delta_prod, x_A, y_A, x_B, y_B = tune_delta(
    x_A, y_A, x_B, y_B,
    L, R_match, Constants, R0_opt, r_grid_prod, psi_prod;
    target_ratio   = 0.5,
    num_tune_steps = num_tune_steps,
    block_size     = 100
)


# Production run
energies_prod, _, _,
E_avg_prod, E_sq_avg_prod, _, _, _, _, _,
_, _ = metropolis(
    N, num_steps_production, delta_prod, L, h, R_match, R0_opt,
    r_grid_prod, psi_prod, u_prime_prod, u_doubleprime_prod, Constants;
    x_init   = vcat(x_A, x_B),
    y_init   = vcat(y_A, y_B),
    progress = true,
    move_all = false
)

n_samples_prod   = length(energies_prod)
global E_production   = E_avg_prod / N
global err_production = sqrt(abs(E_sq_avg_prod - E_avg_prod^2) / n_samples_prod) / N

# Tail correction: AA+BB intra-layer + AB inter-layer
n_layer         = (N/2) / L^2
rc              = L / 2
tail_AA_BB      = 2 * 2π * n_layer / L          # two layers
tail_AB         = π * n_layer * rc^2 / (rc^2 + h^2)^(3/2)
tail_correction = tail_AA_BB + tail_AB
E_corrected     = E_production + tail_correction

println("="^70)
println("PRODUCTION RESULTS  (h = $h r₀)")
println("E/N (raw)       = $(round(E_production,  digits=6)) ± $(round(err_production, digits=6)) ħ²/mr₀²")
println("Tail correction = $(round(tail_correction, digits=6)) ħ²/mr₀²")
println("E/N (corrected) = $(round(E_corrected,   digits=6)) ħ²/mr₀²")
println("="^70 * "\n")

# ──────────────────────────────────────────────────────────────────
# SAVE PRODUCTION DATA
# ──────────────────────────────────────────────────────────────────

prod_path = joinpath(@__DIR__, "..", "data", "production",
                     "VMC_Stage2_N$(N)_nr0sq$(nr0sq)_h$(h).txt")
mkpath(dirname(prod_path))

open(prod_path, "w") do io
    println(io, "# Stage 2: production VMC")
    println(io, "# N=$(N), nr0sq=$(nr0sq), h=$(h), R_match=$(R_match), R0_opt=$(R0_opt)")
    println(io, "# energy_b=$(energy_b_opt)")
    println(io, "# tail_correction=$(tail_correction)")
    println(io, "# E/N (raw)       = $(E_production) ± $(err_production)")
    println(io, "# E/N (corrected) = $(E_corrected)")
    println(io, "#")
    println(io, "# step\tE_local")
    for (i, e) in enumerate(energies_prod)
        println(io, "$(i)\t$(e)")
    end
end

println("✓ Production data saved to: $prod_path")

(; E_production, err_production, tail_correction, E_corrected, prod_path)