# scripts/run_dmc.jl
using DelimitedFiles, Plots, LaTeXStrings, Statistics, ProgressMeter, Base.Threads

include(normpath(joinpath(@__DIR__, "..", "src", "utils.jl")))
include(normpath(joinpath(@__DIR__, "..", "src", "jastrow.jl")))
include(normpath(joinpath(@__DIR__, "..", "src", "energy.jl")))
include(normpath(joinpath(@__DIR__, "..", "src", "metropolis.jl")))
include(normpath(joinpath(@__DIR__, "..", "src", "dmc.jl")))

num_part      = 30
nr0_sq        = 16.0
L             = sqrt(num_part / nr0_sq)
num_steps_dmc = 10^5
num_walkers   = 200
num_target    = 200
Δτ_vals            = [1e-3, 1e-4, 1e-5]

# Read VMC energy as initial E_ref
vmc_path = joinpath(@__DIR__, "..", "data", "results", "VMC",
           "vmc_N$(num_part)_nr0sq$(nr0_sq).txt")
vmc_data = readdlm(vmc_path, '\t', String, skipstart=1)

R_opt = parse(Float64, vmc_data[1, 4])
E_ref_initial = parse(Float64, vmc_data[1, 5])
Error_E_ref_initial = parse(Float64, vmc_data[1, 6])

println("E_ref_initial ± Error_E_ref_initial = ", E_ref_initial, " ± ", Error_E_ref_initial)
println("R_opt  = ", R_opt)

# Get equilibrated VMC config
Constants = calculate_constants(L, R_opt)
x_coord, y_coord = random_initial_config(num_part, L, "Uniform")
delta, x_init, y_init = tune_delta(x_coord, y_coord, L, R_opt, Constants)

# # Load VMC final configuration as initial config for DMC
# config_path = joinpath(@__DIR__, "..", "data", "results",
#               "vmc_config_N$(num_part)_nr0sq$(nr0_sq).txt")
# config_data = readdlm(config_path, '\t', Float64, skipstart=1)
# x_init = config_data[:, 1]
# y_init = config_data[:, 2]
# println("Loaded VMC final config from file")

# Store E_dmc results for different Δτ values
E_dmc_results = Float64[]
E_dmc_err_results = Float64[]
# Run DMC
for Δτ in Δτ_vals
    println("\nRunning DMC with Δτ = ", Δτ)
    E_dmc, E_dmc_err, E_history = dmc(x_init, y_init,
                                     num_walkers, num_part, num_steps_dmc,
                                     Δτ, L, R_opt, Constants,
                                     E_ref_initial, num_target,
                                     plot_energy=true)

    println("\nDMC result for Δτ = $Δτ: E = $E_dmc ± $E_dmc_err")

    push!(E_dmc_results, E_dmc)
    push!(E_dmc_err_results, E_dmc_err)
end


# #Blocking analysis for DMC energy
# block_sizes_dmc = [10, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000]
# sigmas_dmc = Float64[]
# for B in block_sizes_dmc
#     _, sigma_dmc = blocking_statistics(E_history, B)
#     push!(sigmas_dmc, sigma_dmc)
# end

# p = plot(block_sizes_dmc, sigmas_dmc,
#         marker=:circle,
#         xlabel="Block size",
#         ylabel="Standard error",
#         title="Blocking analysis, R_opt = $(round(R_opt, digits=4))",
#         linewidth=2,
#         xticks=block_sizes_dmc,
#         xrotation=45,
#         label="Standard")

# display(p)

# println("Enter plateau block size for standard estimator: ")
# plateau_std = parse(Int, readline())

# avg_E, error_E = blocking_statistics(E_history, plateau_std)
# println("\nDMC energy with blocking analysis: E = $avg_E ± $error_E")
# println("VMC energy: E = $E_ref_initial ± $Error_E_ref_initial")

# Save the results for different values of Δτ
results_path = joinpath(@__DIR__, "..", "data", "results", "DMC",
           "dmc_N$(num_part)_nr0sq$(nr0_sq).txt")
open(results_path, "w") do io
    println(io, "Δτ\tE_dmc\tError_E_dmc")
    for i in 1:length(Δτ_vals)
        println(io, "$(Δτ_vals[i])\t$(E_dmc_results[i])\t$(E_dmc_err_results[i])")
    end

end

# nr0_sq_32 = num_part * nr0_sq^(3/2)
# E_tail = 2π / sqrt(num_part)
# println("E_dmc/N/(nr0^2)^(3/2)          = ", E_dmc/nr0_sq_32, " ± ", E_dmc_err/nr0_sq_32)
# println("E_dmc/N/(nr0^2)^(3/2) + E_tail = ", E_dmc/nr0_sq_32 + E_tail, " ± ", E_dmc_err/nr0_sq_32)