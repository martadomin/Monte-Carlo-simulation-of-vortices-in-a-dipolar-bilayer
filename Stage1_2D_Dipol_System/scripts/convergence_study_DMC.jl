# scripts/run_dmc.jl
using DelimitedFiles, Plots, LaTeXStrings, Statistics, ProgressMeter, Base.Threads

include(normpath(joinpath(@__DIR__, "..", "src", "utils.jl")))
include(normpath(joinpath(@__DIR__, "..", "src", "jastrow.jl")))
include(normpath(joinpath(@__DIR__, "..", "src", "energy.jl")))
include(normpath(joinpath(@__DIR__, "..", "src", "metropolis.jl")))
include(normpath(joinpath(@__DIR__, "..", "src", "dmc.jl")))
include(normpath(joinpath(@__DIR__, "..", "src", "observables.jl")))

# Parameters
num_part      = 30
nr0_sq        = 16.0
L             = sqrt(num_part / nr0_sq)
quadratic     = true

if !quadratic
    type_dmc = "linear"
else
    type_dmc = "quadratic"
end

# Read VMC energy as initial E_ref
vmc_path = joinpath(@__DIR__, "..", "data", "results", "VMC",
           "vmc_N$(num_part)_nr0sq$(nr0_sq).txt")
vmc_data = readdlm(vmc_path, '\t', String, skipstart=1)

R_opt = parse(Float64, vmc_data[1, 4])
E_ref_initial = parse(Float64, vmc_data[1, 5])
Error_E_ref_initial = parse(Float64, vmc_data[1, 6])

println("E_ref_initial ± Error_E_ref_initial = ", E_ref_initial, " ± ", Error_E_ref_initial)
println("R_opt  = ", R_opt)

# Calculate Constants for the Jastrow factor/energy functions
Constants = calculate_constants(L, R_opt)

# Load VMC final configuration as initial config for DMC
config_path = joinpath(@__DIR__, "..", "data", "results", "VMC",
              "vmc_config_N$(num_part)_nr0sq$(nr0_sq).txt")
config_data = readdlm(config_path, '\t', Float64, skipstart=1)
x_init = config_data[:, 1]
y_init = config_data[:, 2]
println("Loaded VMC final config from file")

## Loop over different number of walkers
num_walkers_vals = [200]
## Loop over different time steps
Δτ_vals = [0.01, 0.009, 0.008, 0.007, 0.006, 0.005, 0.004, 0.003, 0.002, 0.001, 0.0009, 0.0008, 0.0007, 0.0006, 0.0005, 0.0004, 0.0003, 0.0002, 0.0001, 0.00009, 0.00008, 0.00007, 0.00006, 0.00005, 0.00004, 0.00003, 0.00002, 0.00001]
num_steps_dmc = 10^5

# Define the path where results will be saved
results_path = joinpath(@__DIR__, "..", "data", "results", "DMC",
           "dmc_N$(num_part)_nr0sq$(nr0_sq)_$(type_dmc).txt")

# Open the file ONCE before the loop begins to write results as they finish
open(results_path, "w") do io
    # Write the header with num_walkers included
    println(io, "num_walkers\tΔτ\tE_dmc\tError_E_dmc")
    
    for num_walkers in num_walkers_vals
        # The target number of walkers is typically the initial number of walkers
        num_target = num_walkers 
        
        for Δτ in Δτ_vals
            println("\nRunning DMC with num_walkers = $num_walkers and Δτ = $Δτ")
            
            # Run DMC, ensuring plot_energy is false to prevent execution pausing
            E_dmc, E_dmc_err, E_history = dmc(x_init, y_init,
                                             num_walkers, num_part, num_steps_dmc,
                                             Δτ, L, R_opt, Constants,
                                             E_ref_initial, num_target,
                                             plot_energy=false,
                                             quadratic=quadratic)

            println("\nRaw DMC result for num_walkers = $num_walkers and Δτ = $Δτ: E = $E_dmc ± $E_dmc_err")

            # Block averaging to get final energy estimates using E_history
            block_sizes = [10, 20, 30, 40, 50, 100, 150, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000]
            sigmas = Float64[]

            for B in block_sizes
                _, sigma = blocking_statistics(E_history, B)
                push!(sigmas, sigma)
            end

            p = plot(block_sizes, sigmas, marker=:circle, label="DMC Error vs Block Size", xlabel="Block Size", ylabel="Error in DMC Energy", title="DMC Error Analysis for num_walkers = $num_walkers and Δτ = $Δτ")
            display(p)

            # --- AUTOMATED PLATEAU DETECTION ---
            println("\n--- Automating Plateau Detection ---")
            plateau_std = detect_plateau(block_sizes, sigmas, window_size=4, rtol=0.05)
            println("Detected plateau block size for standard estimator: ", plateau_std)
            # -----------------------------------

            # Get final error at chosen block size
            avg_energy, sigma = blocking_statistics(E_history, plateau_std)

            println("\nFinal DMC result for num_walkers = $num_walkers and Δτ = $Δτ: E = $avg_energy ± $sigma (using block size = $plateau_std)")

            # Write this specific run's result directly to the text file
            println(io, "$(num_walkers)\t$(Δτ)\t$(avg_energy)\t$(sigma)")
            
            # Flush ensures the line is saved to the hard drive immediately
            flush(io) 
        end
    end
end

println("\nAll DMC runs complete and saved successfully to $results_path")