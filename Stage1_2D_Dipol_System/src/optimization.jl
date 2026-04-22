include(normpath(joinpath(@__DIR__, "jastrow.jl")))
include(normpath(joinpath(@__DIR__, "metropolis.jl")))
include(normpath(joinpath(@__DIR__, "energy.jl")))
include(normpath(joinpath(@__DIR__, "utils.jl")))

using Base.Threads

function sweep_Rmatch(L::Float64, num_part::Int, nr0_sq::Float64, R_match_vals::Vector{Float64},
                      num_steps::Int)::NamedTuple
    energies = Vector{Float64}(undef, length(R_match_vals))
    
    @threads for idx in eachindex(R_match_vals)
        R_match = R_match_vals[idx]
        Constants = calculate_constants(L, R_match)
        x_coord, y_coord = random_initial_config(num_part, L, "Uniform")
        target_ratio = 0.6
        num_tune_steps = 10^4
        block_size = 100
        delta, x_init, y_init = tune_delta(x_coord, y_coord, L, R_match, Constants; target_ratio,
                                            num_tune_steps, block_size)
        E_tot, _, _, _, _ = metropolis(num_part, num_steps, delta, L, R_match, Constants;
                                x_init=x_init, y_init=y_init,
                                final_energy_plot=false,
                                plot_every=10^3)
        energies[idx] = E_tot/(num_part*nr0_sq^(3/2))
    end

    return (R_match_vals=R_match_vals, energies=energies)
end

## Example usage:
# num_part = 30
# nr0_sq = 16.0
# L = sqrt(num_part / nr0_sq)
# R_match_vals = collect(LinRange(0.1*L/2, 0.9*L/2, 2))
# num_steps = 10^6

# _, energies = sweep_Rmatch(L, num_part, nr0_sq, R_match_vals, num_steps)

# print("R_match values: ", R_match_vals, "\n")
# print("Energies per particle: ", energies, "\n")