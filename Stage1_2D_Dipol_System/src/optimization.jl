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
                                                                                                                                
        # # Block averaging to get final energy estimates
        # block_sizes = [10, 20, 30, 40, 50, 100, 150, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1500, 1600, 1700, 1800, 1900, 2000]
        # sigmas = Float64[]
        # sigmas_drift = Float64[]
        # sigmas_laplacian = Float64[]
        # for B in block_sizes
        #     _, sigma = blocking_statistics(energies_vmc, B)
        #     _, sigma_drift = blocking_statistics(energies_drift_vmc, B)
        #     _, sigma_laplacian = blocking_statistics(energies_laplacian_vmc, B)
        #     push!(sigmas, sigma)
        #     push!(sigmas_drift, sigma_drift)
        #     push!(sigmas_laplacian, sigma_laplacian)
        # end

        # # Show plot and ask for plateau block size
        # p = plot(block_sizes, sigmas,
        #          marker=:circle,
        #          xlabel="Block size",
        #          ylabel="Standard error",
        #          title="Blocking analysis, R_match = $(round(R_match, digits=3))",
        #          linewidth=2,
        #          xticks = block_sizes,
        #          xrotation=45)
        # plot!(block_sizes, sigmas_drift,
        #       marker=:square,
        #       linewidth=2)
        # plot!(block_sizes, sigmas_laplacian,
        #       marker=:diamond, 
        #         linewidth=2,
        #         label="Laplacian")
        # display(p)

        # println("\nR_match = $(round(R_match, digits=3))")
        # println("Enter plateau block size for standard estimator: ")
        # plateau_std = parse(Int, readline())
        # println("Enter plateau block size for drift estimator: ")
        # plateau_drift = parse(Int, readline())
        # println("Enter plateau block size for laplacian estimator: ")
        # plateau_laplacian = parse(Int, readline())
        # # Get error at chosen block size
        # avg_energy, sigma = blocking_statistics(energies_vmc, plateau_std)
        # avg_energy_drift, sigma_drift = blocking_statistics(energies_drift_vmc, plateau_drift)
        # avg_energy_laplacian, sigma_laplacian = blocking_statistics(energies_laplacian_vmc, plateau_laplacian)

        # Apply assumed block sizes for error estimation
        # B_std = 400, B_drift = 800, B_laplacian = 1000 are chosen based on previous runs and may be adjusted as needed using the code above for manual selection
        _, sigma = blocking_statistics(energies_vmc, 400)
        _, sigma_drift = blocking_statistics(energies_drift_vmc, 1000)
        _, sigma_laplacian = blocking_statistics(energies_laplacian_vmc, 800)

        # Store averages (normalized by density factor if needed by your plotting script)
        energies[idx] = E_tot
        energies_drift[idx] = E_drift
        energies_laplacian[idx] = E_laplacian

        # Store errors (normalized by density factor)
        error[idx] = sigma
        error_drift[idx] = sigma_drift
        error_laplacian[idx] = sigma_laplacian
    end

    return (R_match_vals=R_match_vals, energies=energies, energies_drift=energies_drift, energies_laplacian=energies_laplacian, error=error, error_drift=error_drift, error_laplacian=error_laplacian)
end