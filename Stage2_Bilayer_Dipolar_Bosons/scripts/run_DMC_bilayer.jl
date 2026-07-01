# scripts/run_dmc.jl

# Read the optimized Rmatch from the Stage I VMC results for the AA and BB contributions.
vmc_path_stage1 = joinpath(@__DIR__, "..", "..", "Stage1_2D_Dipol_System", "data", "results", "VMC",
           "vmc_N$(num_part)_nr0_sq$(nr0_sq / 2).txt")
vmc_data_stage1 = readdlm(vmc_path_stage1, '\t', String, skipstart=1)

R_match = parse(Float64, vmc_data_stage1[1, 4])

# Read the R0 and energy from the Stage II VMC results
vmc_path_stage2 = joinpath(@__DIR__, "..", "..", "Stage2_Bilayer_Dipolar_Bosons", "data", "results", "VMC",
           "VMC_Stage2_N$(num_part)_nr0_sq$(nr0_sq)_h$(h).txt")
vmc_data_stage2 = readdlm(vmc_path_stage2, '\t', String, skipstart=1)

R0 = parse(Float64, vmc_data_stage2[1, 4])
E_ref_initial = parse(Float64, vmc_data_stage2[1, 5]) * num_part

# Calculate Constants for the Jastrow factor/energy functions
Constants = calculate_constants(L, R_match)

# Load VMC final configuration as initial config for DMC
config_path = joinpath(@__DIR__, "..", "data", "configs",
                       "VMC_final_config_N$(num_part)_nr0_sq$(nr0_sq)_h$(h).txt")

config_data = readdlm(config_path, '\t', skipstart=3)  # mixed types -> no Float64 here
layers = String.(config_data[:, 1])
x_init = Float64.(config_data[:, 2])
y_init = Float64.(config_data[:, 3])

# sanity check the A/B split matches your assumed ordering
@assert all(layers[1:num_part÷2] .== "A") "Layer ordering mismatch: expected A particles first"
@assert all(layers[num_part÷2+1:end] .== "B") "Layer ordering mismatch: expected B particles second"

xA_init = x_init[1:num_part÷2]
yA_init = y_init[1:num_part÷2]
xB_init = x_init[num_part÷2+1:end]
yB_init = y_init[num_part÷2+1:end]

println("Loaded VMC final config from file ($(count(==("A"), layers)) A, $(count(==("B"), layers)) B)")

# Define the path where results will be saved
results_path = joinpath(@__DIR__, "..", "..", "Stage2_Bilayer_Dipolar_Bosons", "data", "DMC",
           "dmc_N$(num_part)_nr0_sq$(nr0_sq)_h$(h)_$(type_dmc).txt")

# Open the file ONCE before the loop begins to write results as they finish
open(results_path, "w") do io
    # Write the header with num_walkers included
    println(io, "nr0_sq\th\t\tnum_walkers\tΔτ\tE_dmc\tError_E_dmc")
    
    for num_walkers in num_walkers_vals
        # The target number of walkers is typically the initial number of walkers
        num_target = num_walkers 
        
        for Δτ in Δτ_vals
            println("\nRunning DMC with num_walkers = $num_walkers and Δτ = $Δτ")
            
            # Run DMC, ensuring plot_energy is false to prevent execution pausing
            E_dmc, E_dmc_err, E_history = dmc(xA_init, yA_init,
                                             xB_init, yB_init,
                                             num_walkers, num_part, num_steps_dmc,
                                             Δτ, L, h,
                                             R_match, Constants,
                                             R0, itp_up, itp_upp,
                                             E_ref_initial, num_target,
                                             quadratic=quadratic,
                                             plot_energy=false)

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
            println(io, "$(nr0_sq)\t$(h)\t$(num_walkers)\t$(Δτ)\t$(avg_energy)\t$(sigma)")
            
            # Flush ensures the line is saved to the hard drive immediately
            flush(io) 
        end
    end
end

println("\nAll DMC runs complete and saved successfully to $results_path")