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
r_grid_prod, psi_prod, itp_u_prod, itp_up_prod, itp_upp_prod, _ =
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
    L, R_match, Constants, R0_opt, itp_u_prod;
    target_ratio   = 0.5,
    num_tune_steps = num_tune_steps,
    block_size     = 100
)


# Production run
energies_prod, _, _,
E_avg_prod, E_sq_avg_prod, _, _, _, _, _,
x_coord_final, y_coord_final = metropolis(
    N, num_steps_production, delta_prod, L, h, R_match, R0_opt,
    itp_u_prod, itp_up_prod, itp_upp_prod, Constants;
    x_A_init = x_A, y_A_init = y_A,
    x_B_init = x_B, y_B_init = y_B,
    progress = true,
    move_all = false
)

block_sizes_prod = filter(b -> div(length(energies_prod), b) >= 2,
                          [10,20,30,40,50,100,150,200,300,400,500,600,700,800,900,1000,1200,1500,2000])
sigmas_prod = [let (_, s) = blocking_statistics(energies_prod, B); isnan(s) ? 0.0 : s end
               for B in block_sizes_prod]
plateau_prod = detect_plateau(block_sizes_prod, sigmas_prod; window_size=4, rtol=0.05)
avg_prod, sigma_prod = blocking_statistics(energies_prod, plateau_prod)

global E_production   = avg_prod / N
global err_production = sigma_prod / N


println("="^70)
println("PRODUCTION RESULTS  (h = $h r₀)")
println("E/N (raw)       = $(round(E_production,  digits=6)) ± $(round(err_production, digits=6)) ħ²/mr₀²")
println("="^70 * "\n")

# ──────────────────────────────────────────────────────────────────
# SAVE PRODUCTION DATA
# ──────────────────────────────────────────────────────────────────

prod_path = joinpath(@__DIR__, "..", "data", "production",
                     "VMC_Stage2_N$(N)_nr0sq$(nr0sq)_h$(h).txt")
mkpath(dirname(prod_path))

open(prod_path, "w") do io
    println(io, "N\tnr0_sq\tL\tR_opt\tE/N\tError")
    println(io, "$(N)\t$(nr0sq)\t$(L)\t$(Float64(R0_opt))\t$(E_production)\t$(err_production)")
end

println("✓ Production data saved to: $prod_path")

# ──────────────────────────────────────────────────────────────────
# SAVE FINAL CONFIGURATION (for DMC initialization)
# ──────────────────────────────────────────────────────────────────
config_path = joinpath(@__DIR__, "..", "data", "configs",
                       "VMC_final_config_N$(N)_nr0sq$(nr0sq)_h$(h).txt")
mkpath(dirname(config_path))

open(config_path, "w") do io
    println(io, "# Final VMC configuration for DMC initialization")
    println(io, "# N=$(N), nr0sq=$(nr0sq), h=$(h), L=$(L), R0_opt=$(R0_opt)")
    println(io, "# layer x y   (particles 1:N÷2 = A, N÷2+1:N = B)")
    for i in 1:N
        layer = i <= N÷2 ? "A" : "B"
        println(io, "$layer\t$(x_coord_final[i])\t$(y_coord_final[i])")
    end
end
println("✓ Final configuration saved to: $config_path")

(; E_production, err_production, prod_path, config_path)