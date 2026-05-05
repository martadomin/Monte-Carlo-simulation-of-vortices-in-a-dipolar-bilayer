include(normpath(joinpath(@__DIR__, "jastrow.jl")))
include(normpath(joinpath(@__DIR__, "metropolis.jl")))
include(normpath(joinpath(@__DIR__, "energy.jl")))
include(normpath(joinpath(@__DIR__, "utils.jl")))
include(normpath(joinpath(@__DIR__, "observables.jl")))

using Base.Threads

function sweep_Rmatch(L::Float64, num_part::Int, nr0_sq::Float64, R_match_vals::Vector{Float64},
                      num_steps::Int)::NamedTuple

    energies = Vector{Float64}(undef, length(R_match_vals))
    energies_drift = Vector{Float64}(undef, length(R_match_vals))
    energies_laplacian = Vector{Float64}(undef, length(R_match_vals))
    error = Vector{Float64}(undef, length(R_match_vals))
    error_drift = Vector{Float64}(undef, length(R_match_vals))
    error_laplacian = Vector{Float64}(undef, length(R_match_vals))
    @threads for idx in eachindex(R_match_vals)
        R_match = R_match_vals[idx]
        Constants = calculate_constants(L, R_match)
        x_coord, y_coord = random_initial_config(num_part, L, "Uniform")
        target_ratio = 0.5
        num_tune_steps = 10^4
        block_size = 100
        delta, x_init, y_init = tune_delta(x_coord, y_coord, L, R_match, Constants; target_ratio,
                                            num_tune_steps, block_size)
        energies_vmc, energies_drift_vmc, energies_laplacian_vmc, E_tot, _, E_drift, E_laplacian,E_kin, E_int, _, _, _, _ = metropolis(num_part,
                                                                                                                                num_steps,
                                                                                                                                num_bins=100,
                                                                                                                                delta, L, R_match,
                                                                                                                                Constants;
                                                                                                                                x_init=x_init, y_init=y_init,
                                                                                                                                final_energy_plot=false,
                                                                                                                                plot_every=10^2,
                                                                                                                                progress=true)
                                                                                                                                
        # Block averaging to get final energy estimates
        block_sizes = [10, 20, 30, 40, 50, 100, 150, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1500, 1600, 1700, 1800, 1900, 2000]
        sigmas = Float64[]
        sigmas_drift = Float64[]
        sigmas_laplacian = Float64[]
        for B in block_sizes
            _, sigma = blocking_statistics(energies_vmc, B)
            _, sigma_drift = blocking_statistics(energies_drift_vmc, B)
            _, sigma_laplacian = blocking_statistics(energies_laplacian_vmc, B)
            push!(sigmas, sigma)
            push!(sigmas_drift, sigma_drift)
            push!(sigmas_laplacian, sigma_laplacian)
        end

        # --- AUTOMATED PLATEAU DETECTION ---
        plateau_std = detect_plateau(block_sizes, sigmas, window_size=4, rtol=0.05)
        plateau_drift = detect_plateau(block_sizes, sigmas_drift, window_size=4, rtol=0.05)
        plateau_laplacian = detect_plateau(block_sizes, sigmas_laplacian, window_size=4, rtol=0.05)
        # -----------------------------------
        # Get error at chosen block size
        avg_energy, sigma = blocking_statistics(energies_vmc, plateau_std)
        avg_energy_drift, sigma_drift = blocking_statistics(energies_drift_vmc, plateau_drift)
        avg_energy_laplacian, sigma_laplacian = blocking_statistics(energies_laplacian_vmc, plateau_laplacian)

        # Store averages (normalized by density factor if needed by your plotting script)
        energies[idx] = avg_energy
        energies_drift[idx] = avg_energy_drift
        energies_laplacian[idx] = avg_energy_laplacian

        # Store errors (normalized by density factor)
        error[idx] = sigma
        error_drift[idx] = sigma_drift
        error_laplacian[idx] = sigma_laplacian
    end

    return (R_match_vals=R_match_vals, energies=energies, energies_drift=energies_drift, energies_laplacian=energies_laplacian, error=error, error_drift=error_drift, error_laplacian=error_laplacian)
end