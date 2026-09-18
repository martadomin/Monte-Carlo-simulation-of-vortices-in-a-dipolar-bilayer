# Stage2_Bilayer_Dipolar_Bosons/src/optimization.jl

"""
    load_Rmatch(stage1_dir, N_half, L) -> Float64

Loads R_match optimized for a single layer of N_half particles at box
length L, as saved by Stage1's optimize_Rmatch.jl. Errors if no matching
saved result exists — Stage1's main_VMC.jl must be run first for this
N_half/L combination; this function never runs an optimization itself.
"""
function load_Rmatch(stage1_dir::String, N_half::Int, L::Float64)::Float64
    path = result_path(stage1_dir, "Rmatch_optimum", (N=N_half, L=L))
    if !isfile(path)
        error("No saved R_match found for N=$N_half, L=$L. Run Stage1's main_VMC.jl " *
              "for this N_half/L first — Stage2 does not optimize R_match itself.")
    end
    run = load_run(path; run=1)
    return run.result.R_opt
end

"""
    sweep_R0(h, R_match, Constants, r_min, Δ_shoot, tol_shoot, L, num_part,
             nr0sq, R0_vals, num_steps, stage_dir) -> NamedTuple

For each R0, builds f_AB via build_fAB (getting itp_u/itp_up/itp_upp and
the dimer binding energy energy_b), runs VMC with a BilayerTrial, and
saves the full raw result. R0 values where build_fAB fails (bad ψ, sign
change, near-zero at R0) are recorded as NaN and excluded from the
returned energies — the caller must filter before taking argmin.
"""
function sweep_R0(h::Float64, R_match::Float64, Constants::Tuple{Float64,Float64,Float64},
                   r_min::Float64, Δ_shoot::Float64, tol_shoot::Float64,
                   L::Float64, num_part::Int, nr0sq::Float64,
                   R0_vals::Vector{Float64}, num_steps::Int, stage_dir::String)::NamedTuple

    N_half = num_part ÷ 2
    energies       = Vector{Float64}(undef, length(R0_vals))
    errors         = Vector{Float64}(undef, length(R0_vals))
    energies_b     = Vector{Float64}(undef, length(R0_vals))
    energies_drift = Vector{Float64}(undef, length(R0_vals))
    errors_drift   = Vector{Float64}(undef, length(R0_vals))
    energies_lap   = Vector{Float64}(undef, length(R0_vals))
    errors_lap     = Vector{Float64}(undef, length(R0_vals))

    @threads for i in eachindex(R0_vals)
        R0 = R0_vals[i]
        fAB_result = build_fAB(h, R0, r_min, Δ_shoot, tol_shoot, nr0sq, num_part)

        if fAB_result === nothing
            energies[i] = errors[i] = energies_b[i] = NaN
            energies_drift[i] = errors_drift[i] = energies_lap[i] = errors_lap[i] = NaN
            println("R0 = $(round(R0, digits=4)) → build_fAB failed, skipping")
            continue
        end

        _, _, itp_u, itp_up, itp_upp, energy_b = fAB_result
        trial = BilayerTrial(; R_match, Constants, R0, itp_u, itp_up, itp_upp, h)

        coords = init_random_config((A=N_half, B=N_half), L)
        delta, coords = tune_delta(trial, coords, L; target_ratio=0.5, num_tune_steps=10^4, block_size=100)

        result = metropolis(trial, (A=N_half, B=N_half), num_steps, delta, L;
                             coords_init=coords, progress=false, num_bins=100)

        block_sizes = filter(b -> div(length(result.energies), b) >= 2,
                              [10,20,30,40,50,100,150,200,300,400,500,600,700,800,900,1000,1200,1500,2000])
        sigmas, sigmas_drift, sigmas_lap = Float64[], Float64[], Float64[]
        for B in block_sizes
            _, s  = blocking_statistics(result.energies, B);           push!(sigmas, isnan(s) ? 0.0 : s)
            _, sd = blocking_statistics(result.energies_drift, B);     push!(sigmas_drift, isnan(sd) ? 0.0 : sd)
            _, sl = blocking_statistics(result.energies_laplacian, B); push!(sigmas_lap, isnan(sl) ? 0.0 : sl)
        end
        p  = detect_plateau(block_sizes, sigmas; window_size=4, rtol=0.05)
        pd = detect_plateau(block_sizes, sigmas_drift; window_size=4, rtol=0.05)
        pl = detect_plateau(block_sizes, sigmas_lap; window_size=4, rtol=0.05)

        E_per_N, err               = blocking_statistics(result.energies, p)
        E_drift_per_N, err_drift   = blocking_statistics(result.energies_drift, pd)
        E_lap_per_N, err_lap       = blocking_statistics(result.energies_laplacian, pl)

        energies[i], errors[i], energies_b[i]           = E_per_N/num_part, err/num_part, energy_b
        energies_drift[i], errors_drift[i]               = E_drift_per_N/num_part, err_drift/num_part
        energies_lap[i], errors_lap[i]                   = E_lap_per_N/num_part, err_lap/num_part

        tail = tail_energy(nr0sq, num_part, h)
        println("R0 = $(round(R0, digits=4)) → E/N = $(round(E_per_N/num_part, digits=6)) ± $(round(err/num_part, digits=6)), " *
                "E/N+tail = $(round(E_per_N/num_part + tail, digits=6))")

        path = result_path(stage_dir, "sweep_R0", (N=num_part, nr0sq=nr0sq, h=h, R0=R0))
        save_run(path, (; num_part, nr0sq, h, R_match, R0, num_steps), result)
    end

    return (R0_vals=R0_vals, energies=energies, errors=errors, energies_b=energies_b,
            energies_drift=energies_drift, errors_drift=errors_drift,
            energies_lap=energies_lap, errors_lap=errors_lap)
end

"""
    optimize_R0(h, R_match, Constants, r_min, Δ_shoot, tol_shoot, L, num_part, nr0sq,
                n_points_sweep, num_steps_coarse, num_steps_fine, stage_dir)
        -> (R0_opt, E_opt, err_opt, energy_b_opt)

Coarse sweep over [0, L/2], then a fine sweep over [max(r_min+1e-4, 0.7*R0_rough),
min(L/2, 1.3*R0_rough)]. NaN points (from failed build_fAB) are excluded
before taking the minimum at each stage.
"""
function optimize_R0(h::Float64, R_match::Float64, Constants::Tuple{Float64,Float64,Float64},
                      r_min::Float64, Δ_shoot::Float64, tol_shoot::Float64,
                      L::Float64, num_part::Int, nr0sq::Float64, n_points_sweep::Int,
                      num_steps_coarse::Int, num_steps_fine::Int, stage_dir::String)

    R0_vals_coarse = collect(LinRange(0.0, L/2, n_points_sweep))
    println("\n--- Coarse R0 sweep (h=$h) ---")
    coarse = sweep_R0(h, R_match, Constants, r_min, Δ_shoot, tol_shoot, L, num_part, nr0sq,
                       R0_vals_coarse, num_steps_coarse, stage_dir)

    valid_coarse = .!isnan.(coarse.energies)
    !any(valid_coarse) && error("All R0 points returned NaN in coarse sweep for h=$h")
    idx_rough = findall(valid_coarse)[argmin(coarse.energies[valid_coarse])]
    R0_opt_rough = coarse.R0_vals[idx_rough]
    println("Coarse optimum: R0 = $(round(R0_opt_rough, digits=4)), E/N = $(round(coarse.energies[idx_rough], digits=6))")

    R0_min_fine = max(r_min + 1e-4, 0.7 * R0_opt_rough)
    R0_max_fine = min(L/2, 1.3 * R0_opt_rough)
    R0_vals_fine = collect(LinRange(R0_min_fine, R0_max_fine, n_points_sweep))
    println("\n--- Fine R0 sweep (h=$h) ---")
    fine = sweep_R0(h, R_match, Constants, r_min, Δ_shoot, tol_shoot, L, num_part, nr0sq,
                     R0_vals_fine, num_steps_fine, stage_dir)

    valid_fine = .!isnan.(fine.energies)
    !any(valid_fine) && error("All R0 points returned NaN in fine sweep for h=$h")
    idx_opt = findall(valid_fine)[argmin(fine.energies[valid_fine])]
    R0_opt, E_opt, err_opt, energy_b_opt = fine.R0_vals[idx_opt], fine.energies[idx_opt], fine.errors[idx_opt], fine.energies_b[idx_opt]

    println("\n--- Optimal R0 (h=$h) ---")
    println("R0_opt = $R0_opt, energy_b = $energy_b_opt, E/N = $E_opt ± $err_opt")

    path = result_path(stage_dir, "R0_optimum", (N=num_part, nr0sq=nr0sq, h=h))
    save_run(path, (; num_part, nr0sq, h, R_match, n_points_sweep, num_steps_coarse, num_steps_fine),
             (; coarse, fine, R0_opt_rough, R0_opt, E_opt, err_opt, energy_b_opt))

    return R0_opt, E_opt, err_opt, energy_b_opt
end