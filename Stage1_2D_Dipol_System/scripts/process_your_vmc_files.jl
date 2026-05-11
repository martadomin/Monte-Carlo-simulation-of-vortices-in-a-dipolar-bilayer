#!/usr/bin/env julia
"""
Process VMC results from your specific file format
Reads: num_part  nr0_sq  L  R_opt  E_tot  Error
Computes: E/N and scaled energy Ẽ = (E/N) * (nr0²)^(-3/2)
"""

using DelimitedFiles, Printf

function scale_energy(E_per_N, nr0_sq)
    """Compute scaled energy: (E/N) * (nr₀²)^(-3/2)"""
    return E_per_N * nr0_sq^(-1.5)
end

# ============================================================================
# CONFIGURE YOUR PATHS
# ============================================================================
data_dir = raw"C:\Users\marta\Monte_Carlo_simulation_of_vortices_in_a_dipolar_layer\Stage1_2D_Dipol_System\data\results\VMC"
densities = [16, 32, 48, 64, 96, 128, 196, 256, 384, 512, 768, 1024]

println("\n" * "="^80)
println("VMC ENERGY SCALING - PROCESSING YOUR FILES")
println("="^80)
println("Reading from: $data_dir")
println("="^80)
@printf("%-6s %-4s %-10s %-12s %-12s %-10s %-10s\n", 
        "nr0²", "N", "R_opt", "E_tot", "E/N", "Ẽ", "δẼ")
println("="^80)

results = []
latex_rows = []

for nr0_sq in densities
    # Try with .0 first
    filename = joinpath(data_dir, "vmc_N30_nr0sq$(nr0_sq).0.txt")
    
    # Try without .0 if not found
    if !isfile(filename)
        filename = joinpath(data_dir, "vmc_N30_nr0sq$(nr0_sq).txt")
    end
    
    if !isfile(filename)
        @printf("%-6d ⚠ File not found\n", nr0_sq)
        continue
    end
    
    try
        # Read file (tab-separated, skip header)
        data = readdlm(filename, '\t', skipstart=1)
        
        # Extract values
        # Columns: num_part  nr0_sq  L  R_opt  E_tot  Error
        N = Int(data[1, 1])
        R_opt = data[1, 4]
        E_tot = data[1, 5]
        error = data[1, 6]
        
        # Compute E/N
        E_per_N = E_tot / N
        dE_per_N = error / N
        
        # Compute scaled energy
        E_scaled = scale_energy(E_per_N, nr0_sq)
        dE_scaled = scale_energy(dE_per_N, nr0_sq)
        
        @printf("%-6d %-4d %-10.4f %-12.2f %-12.4f %-10.4f %-10.4f\n",
                nr0_sq, N, R_opt, E_tot, E_per_N, E_scaled, dE_scaled)
        
        push!(results, (nr0_sq, R_opt, E_scaled, dE_scaled))
        push!(latex_rows, @sprintf("%d & %.4f & %.4f & %.4f \\\\", 
                                   nr0_sq, R_opt, E_scaled, dE_scaled))
        
    catch e
        @printf("%-6d ⚠ Error reading file: %s\n", nr0_sq, e)
        continue
    end
end

println("="^80)

if !isempty(latex_rows)
    println("\nLaTeX table rows (ready to paste):")
    println("="^80)
    for row in latex_rows
        println(row)
        println("\\hline")
    end
    println("="^80)
    println("\n✓ Successfully processed $(length(results)) files!")
else
    println("\n⚠ No files were processed. Check your data directory path.")
end
