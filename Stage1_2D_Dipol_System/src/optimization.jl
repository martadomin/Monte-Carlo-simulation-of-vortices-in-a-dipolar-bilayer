function sweep_Rmatch(L::Float64, num_part::Int, R_match_vals::Vector{Float64}, num_steps::Int, stage_dir::String)::NamedTuple

    energies = Vector{Float64}(undef, length(R_match_vals))
    energies_drift = Vector{Float64}(undef, length(R_match_vals))
    energies_laplacian = Vector{Float64}(undef, length(R_match_vals))
    error = Vector{Float64}(undef, length(R_match_vals))
    error_drift = Vector{Float64}(undef, length(R_match_vals))
    error_laplacian = Vector{Float64}(undef, length(R_match_vals))

    @threads for idx in eachindex(R_match_vals)
        R_match = R_match_vals[idx]
        Constants = calculate_constants(L, R_match)
        trial = SingleLayerTrial(; R_match, Constants)

        coords = init_random_config((A=num_part,), L)
        delta, coords = tune_delta(trial, coords, L; target_ratio=0.5, num_tune_steps=10^4, block_size=100)

        result = metropolis(trial, (A=num_part,), num_steps, delta, L;
                             coords_init=coords, final_energy_plot=false,
                             plot_every=10^2, progress=true, num_bins=100)

        block_sizes = [10,20,30,40,50,100,150,200,300,400,500,600,700,800,900,
                       1000,1100,1200,1500,1600,1700,1800,1900,2000]
        sigmas, sigmas_drift, sigmas_laplacian = Float64[], Float64[], Float64[]
        for B in block_sizes
            _, s   = blocking_statistics(result.energies, B)
            _, sd  = blocking_statistics(result.energies_drift, B)
            _, sl  = blocking_statistics(result.energies_laplacian, B)
            push!(sigmas, s); push!(sigmas_drift, sd); push!(sigmas_laplacian, sl)
        end

        plateau_std       = detect_plateau(block_sizes, sigmas, window_size=4, rtol=0.05)
        plateau_drift     = detect_plateau(block_sizes, sigmas_drift, window_size=4, rtol=0.05)
        plateau_laplacian = detect_plateau(block_sizes, sigmas_laplacian, window_size=4, rtol=0.05)

        avg_energy, sigma                   = blocking_statistics(result.energies, plateau_std)
        avg_energy_drift, sigma_drift       = blocking_statistics(result.energies_drift, plateau_drift)
        avg_energy_laplacian, sigma_laplacian = blocking_statistics(result.energies_laplacian, plateau_laplacian)

        energies[idx], error[idx]                     = avg_energy, sigma
        energies_drift[idx], error_drift[idx]         = avg_energy_drift, sigma_drift
        energies_laplacian[idx], error_laplacian[idx] = avg_energy_laplacian, sigma_laplacian

        # Save the FULL raw result for this R_match — not just the summary
        # numbers above. Each R_match gets its own file (distinct path per
        # thread), so this is safe to call from inside @threads with no
        # race condition — no two threads ever write the same file.
        path = result_path(stage_dir, "sweep_Rmatch", (N=num_part, R_match=R_match))
        save_run(path, (; L, num_part, R_match, num_steps), result)
    end

    return (R_match_vals=R_match_vals, energies=energies, energies_drift=energies_drift,
            energies_laplacian=energies_laplacian, error=error,
            error_drift=error_drift, error_laplacian=error_laplacian)
end

"""
    optimize_Rmatch(L, num_part, num_steps_coarse, num_steps_fine, stage_dir) -> Float64

Two-stage R_match optimization: a coarse sweep over the full range, then a
fine sweep centered on the coarse result. Returns R_opt directly — no
intermediate file needs to be written and re-read to get this value.

Saves the combined coarse+fine sweep summary via save_run for later
inspection/plotting, under kind="Rmatch_optimum". Each individual R_match
point's full raw run is already saved separately by sweep_Rmatch itself.
"""
function optimize_Rmatch(L::Float64, num_part::Int, num_steps_coarse::Int,
                          num_steps_fine::Int, stage_dir::String)::Float64

    R_match_vals_coarse = collect(LinRange(0.1*L/2, 0.9*L/2, 16))
    println("\n--- Coarse sweep ($(length(R_match_vals_coarse)) points, $num_steps_coarse steps) ---")
    results_coarse = sweep_Rmatch(L, num_part, R_match_vals_coarse, num_steps_coarse, stage_dir)
    R_opt_rough = results_coarse.R_match_vals[argmin(results_coarse.energies)]
    println("Rough optimal R_match = $R_opt_rough")

    R_match_vals_fine = collect(LinRange(0.7*R_opt_rough, 1.3*R_opt_rough, 16))
    println("\n--- Fine sweep ($(length(R_match_vals_fine)) points, $num_steps_fine steps) ---")
    results_fine = sweep_Rmatch(L, num_part, R_match_vals_fine, num_steps_fine, stage_dir)
    R_opt = results_fine.R_match_vals[argmin(results_fine.energies)]

    println("\n--- Results ---")
    println("Optimal R_match = $R_opt")
    println("Optimal energy per particle = ", minimum(results_fine.energies)/num_part)

    E_opt_total = minimum(results_fine.energies)   # results_fine IS in scope here — it's the function's own local variable

    path = result_path(stage_dir, "Rmatch_optimum", (N=num_part, L=L))
    save_run(path, (; L, num_part, num_steps_coarse, num_steps_fine),
            (; results_coarse, results_fine, R_opt_rough, R_opt, E_opt_total))

    return R_opt
end

"""
    tune_delta_dmc(trial, coords, L; D=0.5, target_ratio=0.5,
                    num_sweep_steps=10^4, max_iters=20) -> (delta_opt, Δτ_DMC, acceptance)

Tunes the Metropolis step size for ALL-PARTICLE moves (move_all=true),
matching DMC's diffusion step where every particle moves simultaneously —
deliberately distinct from tune_delta, which calibrates single-particle
moves for standard VMC. Adjusts `delta` multiplicatively toward
`target_ratio` acceptance, same style as tune_delta, rather than sweeping
a fixed hardcoded list — so it adapts automatically across different
densities/box sizes instead of needing a hand-tuned range per nr0_sq.

Converts the tuned step size to a DMC time step via δ = √(DΔτ) → Δτ = δ²/D.
"""
function tune_delta_dmc(trial::TrialWavefunction, coords::NamedTuple, L::Float64;
                         D::Float64=0.5, target_ratio::Float64=0.5,
                         num_sweep_steps::Int=10^4, max_iters::Int=20)

    delta = 0.1 * L
    acceptance = 0.0

    for iter in 1:max_iters
        result = metropolis(trial, (; (species => length(coords[species].x) for species in species_names(coords))...),
                             num_sweep_steps, delta, L;
                             coords_init=coords, progress=false, move_all=true)
        acceptance = result.acceptance_ratio
        println("Iter $iter: δ = $(round(delta, sigdigits=4)) → acceptance = $(round(acceptance*100, digits=2))%")

        abs(acceptance - target_ratio) < 0.02 && break
        delta *= (acceptance + 1e-2) / target_ratio
        delta = clamp(delta, 1e-6, L/2)
    end

    delta_opt = delta
    Δτ_DMC = delta_opt^2 / (2 * D)
    return delta_opt, Δτ_DMC, acceptance
end