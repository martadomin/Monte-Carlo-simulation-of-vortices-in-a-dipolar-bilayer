"""
optimize_R0.jl

Sweeps R0 values and finds the energy-minimizing interlayer Jastrow cutoff.
All parameters are read from globals set by main_VMC_Stage2.jl:
  N, nr0sq, L, h, R_match, Constants,
  r_min, Δ_shoot, tol_shoot,
  num_steps_coarse, num_steps_fine, num_tune_steps, n_points_sweep

Sets on exit:
  R0_opt, E_opt, err_opt, energy_b_opt
"""

# ──────────────────────────────────────────────────────────────────
# HELPER: run one VMC point for a given R0
# ──────────────────────────────────────────────────────────────────

function run_vmc_point(R0::Float64, num_steps::Int)
    result = build_fAB(h, R0, r_min, Δ_shoot, tol_shoot, nr0sq, N)

    if result === nothing
        @warn "Unphysical solution at R0 = $R0 — skipping"
        return NaN, NaN, NaN
    end

    r_grid, psi, u_prime, u_doubleprime, energy_b = result

    if abs(psi[end] - 1.0) > 1e-6
        @warn "fAB(R0) ≠ 1 — normalization issue at R0 = $R0"
    end

    x_coord, y_coord = random_initial_config(N, L, "Uniform")
    x_A = x_coord[1:N÷2];      y_A = y_coord[1:N÷2]
    x_B = x_coord[N÷2+1:end];  y_B = y_coord[N÷2+1:end]

    delta, x_A, y_A, x_B, y_B = tune_delta(
        x_A, y_A, x_B, y_B,
        L, R_match, Constants, R0, r_grid, psi;
        target_ratio   = 0.5,
        num_tune_steps = num_tune_steps,
        block_size     = 100
    )
    
    energies, _, _,
    E_avg, E_sq_avg, _, _, _, _, _,
    _, _ = metropolis(
        N, num_steps, delta, L, h, R_match, R0, r_grid, psi,
        u_prime, u_doubleprime, Constants;
        x_init   = vcat(x_A, x_B),
        y_init   = vcat(y_A, y_B),
        progress = false,
        move_all = false
    )

    n_samples = length(energies)
    E_per_N   = E_avg / N
    error     = sqrt(abs(E_sq_avg - E_avg^2) / n_samples) / N

    return E_per_N, error, energy_b
end

# ──────────────────────────────────────────────────────────────────
# COARSE SWEEP
# ──────────────────────────────────────────────────────────────────

R0_min_coarse = 0.5 * h
R0_max_coarse = L/2
R0_vals_coarse = collect(LinRange(R0_min_coarse, R0_max_coarse, n_points_sweep))

energies_coarse = Vector{Float64}(undef, n_points_sweep)
errors_coarse   = Vector{Float64}(undef, n_points_sweep)
eps_b_coarse    = Vector{Float64}(undef, n_points_sweep)

println("="^70)
println("COARSE SWEEP  (h = $h r₀)")
println("R0 range : [$(round(R0_min_coarse,digits=4)), $(round(R0_max_coarse,digits=4))] r₀")
println("L/2      : $(round(L/2, digits=4)) r₀")
println("N points : $n_points_sweep")
println("MC steps : $num_steps_coarse")
println("Threads  : $(Threads.nthreads())")
println("="^70 * "\n")

@threads for i in eachindex(R0_vals_coarse)
    R0 = R0_vals_coarse[i]
    println("[$i/$n_points_sweep] R0 = $(round(R0, digits=4)) r₀  (thread $(Threads.threadid()))")
    E_per_N, error, energy_b = run_vmc_point(R0, num_steps_coarse)
    energies_coarse[i] = E_per_N
    errors_coarse[i]   = error
    eps_b_coarse[i]    = energy_b
    println("E/N = $(round(E_per_N, digits=6)) ± $(round(error, digits=6)) ħ²/mr₀²\n")
end

_, idx_rough = findmin(energies_coarse)
R0_opt_rough = R0_vals_coarse[idx_rough]
E_opt_rough  = energies_coarse[idx_rough]

println("="^70)
println("Coarse optimum: R0 = $(round(R0_opt_rough, digits=4)) r₀,  E/N = $(round(E_opt_rough, digits=6))")
println("="^70 * "\n")

# ──────────────────────────────────────────────────────────────────
# FINE SWEEP
# ──────────────────────────────────────────────────────────────────

R0_min_fine  = max(r_min + 1e-4, 0.7 * R0_opt_rough)
R0_max_fine  = min(L/2,          1.3 * R0_opt_rough)
R0_vals_fine = collect(LinRange(R0_min_fine, R0_max_fine, n_points_sweep))

energies_fine = Vector{Float64}(undef, n_points_sweep)
errors_fine   = Vector{Float64}(undef, n_points_sweep)
eps_b_fine    = Vector{Float64}(undef, n_points_sweep)

println("="^70)
println("FINE SWEEP  (h = $h r₀)")
println("R0 range : [$(round(R0_min_fine,digits=4)), $(round(R0_max_fine,digits=4))] r₀")
println("N points : $n_points_sweep")
println("MC steps : $num_steps_fine")
println("Threads  : $(Threads.nthreads())")
println("="^70 * "\n")

@threads for i in eachindex(R0_vals_fine)
    R0 = R0_vals_fine[i]
    println("[$i/$n_points_sweep] R0 = $(round(R0, digits=4)) r₀  (thread $(Threads.threadid()))")
    E_per_N, error, energy_b = run_vmc_point(R0, num_steps_fine)
    energies_fine[i] = E_per_N
    errors_fine[i]   = error
    eps_b_fine[i]    = energy_b
    println("    E/N = $(round(E_per_N, digits=6)) ± $(round(error, digits=6)) ħ²/mr₀²\n")
end

_, idx_opt = findmin(energies_fine)
global R0_opt       = R0_vals_fine[idx_opt]
global E_opt        = energies_fine[idx_opt]
global err_opt      = errors_fine[idx_opt]
global energy_b_opt = eps_b_fine[idx_opt]

println("="^70)
println("OPTIMAL R0  (h = $h r₀)")
println("R0_opt = $(round(R0_opt,       digits=6)) r₀")
println("ε_b    = $(round(energy_b_opt, digits=6))")
println("E/N    = $(round(E_opt,        digits=6)) ± $(round(err_opt, digits=6)) ħ²/mr₀²")
println("="^70 * "\n")

# ──────────────────────────────────────────────────────────────────
# SAVE SWEEP RESULTS
# ──────────────────────────────────────────────────────────────────

sweep_path = joinpath(@__DIR__, "..", "data", "sweep_results",
                      "R0_sweep_Stage2_N$(N)_nr0sq$(nr0sq)_h$(h).txt")
mkpath(dirname(sweep_path))

open(sweep_path, "w") do io
    println(io, "# Stage 2: R0 optimization sweep")
    println(io, "# N=$(N), nr0sq=$(nr0sq), h=$(h), R_match=$(R_match)")
    println(io, "# L=$(L), L/2=$(L/2)")
    println(io, "#")
    println(io, "# --- Coarse sweep ---")
    println(io, "# R0\tE/N\tError\teps_b")
    for (r0, e, err, eb) in zip(R0_vals_coarse, energies_coarse, errors_coarse, eps_b_coarse)
        println(io, "$(r0)\t$(e)\t$(err)\t$(eb)")
    end
    println(io, "#")
    println(io, "# --- Fine sweep ---")
    println(io, "# R0\tE/N\tError\teps_b")
    for (r0, e, err, eb) in zip(R0_vals_fine, energies_fine, errors_fine, eps_b_fine)
        println(io, "$(r0)\t$(e)\t$(err)\t$(eb)")
    end
    println(io, "#")
    println(io, "# --- Optimal ---")
    println(io, "R0_opt      = $(R0_opt)")
    println(io, "energy_b    = $(energy_b_opt)")
    println(io, "E_opt       = $(E_opt)")
    println(io, "err_opt     = $(err_opt)")
end

println("✓ Sweep results saved to: $sweep_path")