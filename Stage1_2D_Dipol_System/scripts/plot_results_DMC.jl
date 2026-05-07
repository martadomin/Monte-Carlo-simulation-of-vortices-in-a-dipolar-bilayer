using DelimitedFiles, Plots, LaTeXStrings

ENV["PATH"] = "C:\\Users\\marta\\AppData\\Local\\Programs\\MiKTeX\\miktex\\bin\\x64;" * ENV["PATH"]
pgfplotsx()

num_part  = 30
nr0_sq    = 16.0
nr0_sq_32 = num_part * nr0_sq^(3/2)
E_tail    = 2π / sqrt(num_part)
quadratic = false

if !quadratic
    type_dmc = "linear"
else
    type_dmc = "quadratic"
end 

# Load data Linear DMC
dmc_path_lin = joinpath(@__DIR__, "..", "data", "results", "DMC",
           "dmc_N$(num_part)_nr0sq$(nr0_sq)_linear.txt")
dmc_data_lin = readdlm(dmc_path_lin, '\t', Float64, skipstart=1)

num_walkers_vals_lin = dmc_data_lin[:, 1]
Δτ_vals_lin          = dmc_data_lin[:, 2]
E_dmc_lin            = dmc_data_lin[:, 3]
Error_E_dmc_lin      = dmc_data_lin[:, 4]

#Load data Quadratic DMC
dmc_path_quad = joinpath(@__DIR__, "..", "data", "results", "DMC",
              "dmc_N$(num_part)_nr0sq$(nr0_sq)_quadratic.txt")
dmc_data_quad = readdlm(dmc_path_quad, '\t', Float64, skipstart=1)

num_walkers_vals_quad = dmc_data_quad[:, 1]
Δτ_vals_quad          = dmc_data_quad[:, 2]
E_dmc_quad            = dmc_data_quad[:, 3]
Error_E_dmc_quad      = dmc_data_quad[:, 4]

#Load data VMC
vmc_path = joinpath(@__DIR__, "..", "data", "results", "VMC",
                    "vmc_N$(num_part)_nr0sq$(nr0_sq).txt")
vmc_data    = readdlm(vmc_path, '\t', Float64, skipstart=1)
E_vmc       = vmc_data[1, 5]
Error_E_vmc = vmc_data[1, 6]

#Load data 

E_norm       = E_dmc ./ nr0_sq_32
Error_norm   = Error_E_dmc ./ nr0_sq_32

unique_walkers = sort(unique(num_walkers_vals))
unique_Δτ      = sort(unique(Δτ_vals))
colors         = palette(:viridis, length(unique_walkers))

# ------------------------------------------
# Plot 1: Energy vs Δτ for each num_walkers
# ------------------------------------------
p1 = plot(
    xlabel  = L"\Delta\tau",
    ylabel  = L"E/N \cdot (nr_0^2)^{-3/2}",
    title   = L"DMC\ \mathrm{Energy\ vs}\ \Delta\tau,\ N=30,\ nr_0^2=16",
    legend  = :topright,
    xscale  = :log10,
    framestyle = :box,
    grid    = true,
    gridalpha = 0.3
)

for (idx, nw) in enumerate(unique_walkers)
    mask = num_walkers_vals .== nw
    Δτ_sub = Δτ_vals[mask]
    E_sub  = E_norm[mask]
    Err_sub = Error_norm[mask]
    
    plot!(p1, Δτ_sub, E_sub,
          yerror   = Err_sub,
          label    = L"N_w = %$(Int(nw))",
          marker   = :circle,
          markersize = 5,
          linewidth  = 2,
          color    = colors[idx])
end

hline!(p1, [(E_vmc - Error_E_vmc) / nr0_sq_32,
            (E_vmc + Error_E_vmc) / nr0_sq_32];
       label     = "",
       linestyle = :dot,
       linewidth = 1,
       color     = :red,
       alpha     = 0.5)

savefig(p1, joinpath(@__DIR__, "..", "data", "plots", "DMC", "dmc_E_vs_dtau_N$(num_part)_nr0sq$(nr0_sq).pdf"))
println("Saved Plot 1")
display(p1)

# # ------------------------------------------
# # Plot 2: Energy vs num_walkers for each Δτ
# # ------------------------------------------
# colors2 = palette(:plasma, length(unique_Δτ))

# p2 = plot(
#     xlabel  = L"N_{\mathrm{walkers}}",
#     ylabel  = L"E/N \cdot (nr_0^2)^{-3/2}",
#     title   = L"DMC\ \mathrm{Energy\ vs}\ N_{\mathrm{walkers}},\ N=30,\ nr_0^2=16",
#     legend  = :topright,
#     framestyle = :box,
#     grid    = true,
#     gridalpha = 0.3
# )

# for (idx, dt) in enumerate(unique_Δτ)
#     mask = Δτ_vals .== dt
#     nw_sub  = num_walkers_vals[mask]
#     E_sub   = E_norm[mask]
#     Err_sub = Error_norm[mask]
    
#     plot!(p2, nw_sub, E_sub,
#           yerror   = Err_sub,
#           label    = L"\Delta\tau = %$(dt)",
#           marker   = :circle,
#           markersize = 5,
#           linewidth  = 2,
#           color    = colors2[idx])
# end

# savefig(p2, joinpath(@__DIR__, "..", "data", "plots", "DMC", "dmc_E_vs_walkers_N$(num_part)_nr0sq$(nr0_sq).pdf"))
# println("Saved Plot 2")
# display(p2)

