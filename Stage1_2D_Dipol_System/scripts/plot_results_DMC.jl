using DelimitedFiles, Plots, LaTeXStrings

ENV["PATH"] = "C:\\Users\\marta\\AppData\\Local\\Programs\\MiKTeX\\miktex\\bin\\x64;" * ENV["PATH"]
gr()

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
           "dmc_N$(num_part)_nr0sq$(nr0_sq)_linear_no_branching.txt")
dmc_data_lin = readdlm(dmc_path_lin, '\t', Float64, skipstart=1)

num_walkers_vals_lin = dmc_data_lin[:, 1]
Δτ_vals_lin          = dmc_data_lin[:, 2]
E_dmc_lin            = dmc_data_lin[:, 3]
Error_E_dmc_lin      = dmc_data_lin[:, 4]

#Load data Quadratic DMC
dmc_path_quad = joinpath(@__DIR__, "..", "data", "results", "DMC",
              "dmc_N$(num_part)_nr0sq$(nr0_sq)_quadratic_no_branching.txt")
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

E_norm_lin      = E_dmc_lin ./ nr0_sq_32
Error_norm_lin   = Error_E_dmc_lin ./ nr0_sq_32

E_norm_quad = E_dmc_quad ./ nr0_sq_32
Error_norm_quad = Error_E_dmc_quad ./ nr0_sq_32


unique_walkers = sort(unique(num_walkers_vals_lin))
unique_Δτ      = sort(unique(Δτ_vals_lin))
colors         = palette(:viridis, length(unique_walkers)+1)

# ------------------------------------------
# Plot 1: Energy vs Δτ for each num_walkers
# ------------------------------------------
p1 = plot(
    xlabel  = L"\Delta\tau",
    ylabel  = L"E/N \cdot (nr_0^2)^{-3/2}",
    title   = L"DMC\ \mathrm{Energy\ vs}\ \Delta\tau,\ N=30,\ nr_0^2=16",
    legend  = :topleft,
    framestyle = :box,
    grid    = true,
    gridalpha = 0.3,
    xlim = (1e-5, 1e-4),
    ylim = (5.64, 5.69)
)

for (idx, nw) in enumerate(unique_walkers)
    # separate masks for each file
    mask_lin  = num_walkers_vals_lin  .== nw
    mask_quad = num_walkers_vals_quad .== nw

    Δτ_sub_lin   = Δτ_vals_lin[mask_lin]
    E_sub_lin    = E_norm_lin[mask_lin]
    Err_sub_lin  = Error_norm_lin[mask_lin]    # ← was Error_norm_quad

    Δτ_sub_quad  = Δτ_vals_quad[mask_quad]
    E_sub_quad   = E_norm_quad[mask_quad]
    Err_sub_quad = Error_norm_quad[mask_quad]

    plot!(p1, Δτ_sub_lin, E_sub_lin;
          yerror = Err_sub_lin,
          label  = L"N_w = %$(Int(nw)) \; \mathrm{linear}",
          marker = :circle, markersize = 5, linewidth = 2,
          color  = colors[idx])

    plot!(p1, Δτ_sub_quad, E_sub_quad;
          yerror    = Err_sub_quad,
          label     = L"N_w = %$(Int(nw)) \; \mathrm{quadratic}",
          marker    = :square, markersize = 5, linewidth = 2,
          color     = colors[idx+1], linestyle = :dash)
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

