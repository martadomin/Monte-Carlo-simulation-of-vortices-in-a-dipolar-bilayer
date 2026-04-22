include("metropolis_VMC/metropolis.jl")

#System parameters
# n = [16, 32, 48, 64, 96, 128, 196, 256] # Density values to sweep over
nr0_sq = 16.0
num_part = 30
L = sqrt(num_part/nr0_sq)

println("Density n = ", nr0_sq)
println("Number of particles: ", num_part)
println("L = ", L)

#Simulation parameters
R_match_vals = LinRange(0.1*L/2, 0.9*L/2, 20) # Rough sweep over R_match values to find optimal one for energy minimization
println("R_match values to sweep over: ", R_match_vals)

results_path = joinpath(@__DIR__, "Data_VMC/R_match_energy_results_N$(num_part)_nr0sq$(nr0_sq).txt")
R_match_vals_vec = collect(R_match_vals)
energy_means = Vector{Float64}(undef, length(R_match_vals_vec))
energy_sq_means = Vector{Float64}(undef, length(R_match_vals_vec))
acceptance_ratio_vec = Vector{Float64}(undef, length(R_match_vals_vec))
energy_kin_means = Vector{Float64}(undef, length(R_match_vals_vec))
energy_int_means = Vector{Float64}(undef, length(R_match_vals_vec))

print("Running VMC sweep over R_match values... ")
#Start timer for the entire sweep
@time begin
    @threads for idx in eachindex(R_match_vals_vec)
        R_match = R_match_vals_vec[idx]
        num_steps = 10^6
        num_bins = 100
        delta = 0.1
        Constants = calculate_constants(L, R_match)
        E_mean, E_sq_mean, acceptance_ratio, E_kin_mean, E_int_mean = metropolis(num_part, num_steps, num_bins, delta, L, R_match, Constants; live_plot=false, plot_every=10^3)
        energy_means[idx] = E_mean/(num_part*nr0_sq^(3/2))
        energy_sq_means[idx] = E_sq_mean
        acceptance_ratio_vec[idx] = acceptance_ratio
        energy_kin_means[idx] = E_kin_mean
        energy_int_means[idx] = E_int_mean
    end

    open(results_path, "w") do io
        println(io, "R_match\tMeanEnergyPerParticle\tMeanSquaredEnergy\tKineticEnergy\tPotentialEnergy\tAcceptanceRatio")
        for idx in eachindex(R_match_vals_vec)
            println(io, "$(R_match_vals_vec[idx])\t$(energy_means[idx])\t$(energy_sq_means[idx])\t$(energy_kin_means[idx])\t$(energy_int_means[idx])\t$(acceptance_ratio_vec[idx])")
        end
    end
end

println("VMC sweep completed.")
println("Saved sweep results to: ", results_path)
