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
        return NaN, NaN, NaN, NaN, NaN, NaN, NaN
    end

    r_grid, psi, itp_u, itp_up, itp_upp, energy_b = result

    x_coord, y_coord = random_initial_config(N, L, "Uniform")
    x_A = x_coord[1:N÷2];      y_A = y_coord[1:N÷2]
    x_B = x_coord[N÷2+1:end];  y_B = y_coord[N÷2+1:end]

    delta, x_A, y_A, x_B, y_B = tune_delta(
        x_A, y_A, x_B, y_B,
        L, R_match, Constants, R0, itp_u;
        target_ratio   = 0.5,
        num_tune_steps = num_tune_steps,
        block_size     = 100
    )

    energies, energies_drift, energies_laplacian,
    E_avg, E_sq_avg, E_avg_drift, E_avg_laplacian, _, _, _,
    _, _ = metropolis(
        N, num_steps, delta, L, h, R_match, R0, itp_u,
        itp_up, itp_upp, Constants;
        x_A_init = x_A, y_A_init = y_A,
        x_B_init = x_B, y_B_init = y_B,
        progress = false,
        move_all = false
    )

    n_samples = length(energies)
    E_per_N         = E_avg / N
    error           = sqrt(abs(E_sq_avg - E_avg^2) / n_samples) / N
    E_drift_per_N   = E_avg_drift / N
    error_drift     = sqrt(abs(sum(abs2, energies_drift)/n_samples - E_avg_drift^2) / n_samples) / N
    E_lap_per_N     = E_avg_laplacian / N
    error_laplacian = sqrt(abs(sum(abs2, energies_laplacian)/n_samples - E_avg_laplacian^2) / n_samples) / N

    return E_per_N, error, energy_b, E_drift_per_N, error_drift, E_lap_per_N, error_laplacian
end

# ──────────────────────────────────────────────────────────────────
# COARSE SWEEP
# ──────────────────────────────────────────────────────────────────

R0_min_coarse = 0.5 * h
R0_max_coarse = min(L/2, 7.0 * h)
R0_vals_coarse = collect(LinRange(R0_min_coarse, R0_max_coarse, n_points_sweep))

energies_coarse   = Vector{Float64}(undef, n_points_sweep)
errors_coarse     = Vector{Float64}(undef, n_points_sweep)
eps_b_coarse      = Vector{Float64}(undef, n_points_sweep)
energies_drift_coarse  = Vector{Float64}(undef, n_points_sweep)
errors_drift_coarse    = Vector{Float64}(undef, n_points_sweep)
energies_lap_coarse    = Vector{Float64}(undef, n_points_sweep)
errors_lap_coarse      = Vector{Float64}(undef, n_points_sweep)

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
    E_per_N, error, energy_b, E_drift, err_drift, E_lap, err_lap = run_vmc_point(R0, num_steps_coarse)
    energies_coarse[i]       = E_per_N
    errors_coarse[i]         = error
    eps_b_coarse[i]          = energy_b
    energies_drift_coarse[i] = E_drift
    errors_drift_coarse[i]   = err_drift
    energies_lap_coarse[i]   = E_lap
    errors_lap_coarse[i]     = err_lap
    println("E/N = $(round(E_per_N, digits=6)) ± $(round(error, digits=6)) ħ²/mr₀²\n")
end

valid = .!isnan.(energies_coarse)
if !any(valid)
    error("All R0 points returned NaN in coarse sweep for h=$h")
end

_, idx_rough_valid = findmin(energies_coarse[valid])
idx_rough    = findall(valid)[idx_rough_valid]
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

energies_fine   = Vector{Float64}(undef, n_points_sweep)
errors_fine     = Vector{Float64}(undef, n_points_sweep)
eps_b_fine      = Vector{Float64}(undef, n_points_sweep)
energies_drift_fine = Vector{Float64}(undef, n_points_sweep)
errors_drift_fine   = Vector{Float64}(undef, n_points_sweep)
energies_lap_fine   = Vector{Float64}(undef, n_points_sweep)
errors_lap_fine     = Vector{Float64}(undef, n_points_sweep)

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
    E_per_N, error, energy_b, E_drift, err_drift, E_lap, err_lap = run_vmc_point(R0, num_steps_fine)
    energies_fine[i]       = E_per_N
    errors_fine[i]         = error
    eps_b_fine[i]          = energy_b
    energies_drift_fine[i] = E_drift
    errors_drift_fine[i]   = err_drift
    energies_lap_fine[i]   = E_lap
    errors_lap_fine[i]     = err_lap
    println("    E/N = $(round(E_per_N, digits=6)) ± $(round(error, digits=6)) ħ²/mr₀²\n")
end

valid_fine = .!isnan.(energies_fine)
if !any(valid_fine)
    error("All R0 points returned NaN in fine sweep for h=$h")
end

_, idx_opt_valid    = findmin(energies_fine[valid_fine])
idx_opt             = findall(valid_fine)[idx_opt_valid]

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
    println(io, "# R0\tE/N\tError\teps_b\tE_drift/N\tError_drift\tE_lap/N\tError_lap")
    for (r0, e, err, eb, ed, eed, el, eel) in zip(
            R0_vals_coarse, energies_coarse, errors_coarse, eps_b_coarse,
            energies_drift_coarse, errors_drift_coarse,
            energies_lap_coarse, errors_lap_coarse)
        println(io, "$(r0)\t$(e)\t$(err)\t$(eb)\t$(ed)\t$(eed)\t$(el)\t$(eel)")
    end
    println(io, "#")
    println(io, "# --- Fine sweep ---")
    println(io, "# R0\tE/N\tError\teps_b\tE_drift/N\tError_drift\tE_lap/N\tError_lap")
    for (r0, e, err, eb, ed, eed, el, eel) in zip(
            R0_vals_fine, energies_fine, errors_fine, eps_b_fine,
            energies_drift_fine, errors_drift_fine,
            energies_lap_fine, errors_lap_fine)
        println(io, "$(r0)\t$(e)\t$(err)\t$(eb)\t$(ed)\t$(eed)\t$(el)\t$(eel)")
    end
    println(io, "#")
    println(io, "# --- Optimal ---")
    println(io, "R0_opt      = $(R0_opt)")
    println(io, "energy_b    = $(energy_b_opt)")
    println(io, "E_opt       = $(E_opt)")
    println(io, "err_opt     = $(err_opt)")
end

println("✓ Sweep results saved to: $sweep_path")