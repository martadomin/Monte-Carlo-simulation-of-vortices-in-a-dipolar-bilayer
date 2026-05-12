using DelimitedFiles, Plots, LaTeXStrings, Polynomials

ENV["PATH"] = "C:\\Users\\marta\\AppData\\Local\\Programs\\MiKTeX\\miktex\\bin\\x64;" * ENV["PATH"]
# pgfplotsx()
gr()

num_part  = 30
nr0_sq    = 16.0
nr0_sq_32 = num_part * nr0_sq^(3/2)
E_tail    = 2π / sqrt(num_part)
quadratic = false
num_walkers = 200

if !quadratic
    type_dmc = "linear"
else
    type_dmc = "quadratic"
end 

# Load data Linear DMC
dmc_path_lin = joinpath(@__DIR__, "..", "data", "results", "DMC",
           "dmc_N$(num_part)_nr0sq$(nr0_sq)_Nw$(num_walkers)_linear_no_branching.txt")
dmc_data_lin = readdlm(dmc_path_lin, '\t', Float64, skipstart=1)

num_walkers_vals_lin = dmc_data_lin[:, 1]
Δτ_vals_lin          = dmc_data_lin[:, 2]
E_dmc_lin            = dmc_data_lin[:, 3]
Error_E_dmc_lin      = dmc_data_lin[:, 4]

# Load data Linear DMC
dmc_path_lin_real = joinpath(@__DIR__, "..", "data", "results", "DMC",
           "dmc_N$(num_part)_nr0sq$(nr0_sq)_Nw$(num_walkers)_linear.txt")
dmc_data_lin_real = readdlm(dmc_path_lin_real, '\t', Float64, skipstart=1)

num_walkers_vals_lin_real = dmc_data_lin_real[:, 1]
Δτ_vals_lin_real          = dmc_data_lin_real[:, 2]
E_dmc_lin_real            = dmc_data_lin_real[:, 3]
Error_E_dmc_lin_real      = dmc_data_lin_real[:, 4]

#Load data Quadratic DMC
dmc_path_quad = joinpath(@__DIR__, "..", "data", "results", "DMC",
              "dmc_N$(num_part)_nr0sq$(nr0_sq)_Nw$(num_walkers)_quadratic_no_branching.txt")
dmc_data_quad = readdlm(dmc_path_quad, '\t', Float64, skipstart=1)

num_walkers_vals_quad = dmc_data_quad[:, 1]
Δτ_vals_quad = dmc_data_quad[:, 2]
E_dmc_quad = dmc_data_quad[:, 3]
Error_E_dmc_quad = dmc_data_quad[:, 4]

#Load data Quadratic DMC
dmc_path_quad_real = joinpath(@__DIR__, "..", "data", "results", "DMC",
              "dmc_N$(num_part)_nr0sq$(nr0_sq)_Nw$(num_walkers)_quadratic.txt")
dmc_data_quad_real = readdlm(dmc_path_quad_real, '\t', Float64, skipstart=1)

num_walkers_vals_quad_real = dmc_data_quad_real[:, 1]
Δτ_vals_quad_real          = dmc_data_quad_real[:, 2]
E_dmc_quad_real            = dmc_data_quad_real[:, 3]
Error_E_dmc_quad_real      = dmc_data_quad_real[:, 4]


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

E_norm_quad_real = E_dmc_quad_real ./ nr0_sq_32
Error_norm_quad_real = Error_E_dmc_quad_real ./ nr0_sq_32

E_norm_lin_real = E_dmc_lin_real ./ nr0_sq_32
Error_norm_lin_real = Error_E_dmc_lin_real ./ nr0_sq_32

unique_walkers = sort(unique(num_walkers_vals_lin_real))
unique_Δτ      = sort(unique(Δτ_vals_lin_real))
colors         = palette(:viridis, length(unique_walkers)+1)

# ------------------------------------------
# Plot 1: Energy vs Δτ for each num_walkers
# ------------------------------------------
p1 = plot(
    xlabel  = L"$\Delta\tau$",
    ylabel  = L"$E/N \cdot (nr_0^2)^{-3/2}$",
    title   = L"$\mathrm{DMC\ Energy\ vs}\ \Delta\tau,\ N=30,\ nr_0^2=16$",
    legend  = :outerright,
    framestyle = :box,
    tickfontfamily="Computer Modern",
    grid    = true,
    dpi = 600
)

for (idx, nw) in enumerate(unique_walkers)
    mask_lin  = (num_walkers_vals_lin  .== nw) .& (Δτ_vals_lin  .>= 1e-5) .& (Δτ_vals_lin  .<= 1e-3)
    mask_quad = (num_walkers_vals_quad .== nw) .& (Δτ_vals_quad .>= 1e-5) .& (Δτ_vals_quad .<= 1e-3)
    mask_lin_real  = (num_walkers_vals_lin_real  .== nw) .& (Δτ_vals_lin_real  .>= 1e-5) .& (Δτ_vals_lin_real  .<= 1e-3)
    mask_quad_real = (num_walkers_vals_quad_real .== nw) .& (Δτ_vals_quad_real .>= 1e-5) .& (Δτ_vals_quad_real .<= 1e-3)

    Δτ_sub_lin  = Δτ_vals_lin[mask_lin]
    E_sub_lin   = E_norm_lin[mask_lin]
    err_sub_lin = Error_norm_lin[mask_lin]

    Δτ_sub_quad = Δτ_vals_quad[mask_quad]
    E_sub_quad  = E_norm_quad[mask_quad]
    err_sub_quad = Error_norm_quad[mask_quad]

    Δτ_sub_lin_real  = Δτ_vals_lin_real[mask_lin_real]
    E_sub_lin_real   = E_norm_lin_real[mask_lin_real]
    err_sub_lin_real = Error_norm_lin_real[mask_lin_real]

    Δτ_sub_quad_real = Δτ_vals_quad_real[mask_quad_real]
    E_sub_quad_real  = E_norm_quad_real[mask_quad_real]
    err_sub_quad_real = Error_norm_quad_real[mask_quad_real]

    # Data points
    plot!(p1, Δτ_sub_lin, E_sub_lin;
      label = L"$N_w = %$nw\ \mathrm{(linear\ no\ branching)}$",
      color = :royalblue, 
      alpha = 0.7, 
      linestyle = :dashdot, # Distinct from the fit line
      marker = :circle, 
      markersize = 5, 
      markerstrokewidth = 0.5)
    
    plot!(p1, Δτ_sub_lin_real, E_sub_lin_real;
      label = L"$N_w = %$nw\ \mathrm{(linear)}$",
      color = :midnightblue, 
      alpha = 0.7, 
      linestyle = :solid, # Distinct from the fit line
      marker = :circle, 
      markersize = 5, 
      markerstrokewidth = 0.5)

    plot!(p1, Δτ_sub_quad, E_sub_quad;
      label = L"$N_w = %$nw\ \mathrm{(quadratic\ no\ branching)}$",
      color = :chocolate, 
      linewidth = 2, 
      marker = :diamond, 
      markersize = 6,
      markerstrokewidth = 0.5,
      linestyle = :dash)
    
    plot!(p1, Δτ_sub_quad_real, E_sub_quad_real;
      label = L"$N_w = %$nw\ \mathrm{(quadratic\ real)}$",
      color = :saddlebrown,
      linewidth = 2,
      marker = :diamond,
      markersize = 6,
      markerstrokewidth = 0.5,
      linestyle = :solid)

    # Linear fit to linear DMC:  E = a₀ + a₁Δτ
    fit_lin  = fit(Δτ_sub_lin,  E_sub_lin,  1)
    fit_lin_real  = fit(Δτ_sub_lin_real,  E_sub_lin_real,  1)

    # Quadratic fit to quadratic DMC:  E = b₀ + b₁Δτ + b₂(Δτ)²
    fit_quad = fit(Δτ_sub_quad, E_sub_quad, 2)
    fit_quad_real = fit(Δτ_sub_quad_real, E_sub_quad_real, 2)

    Δτ_range = range(0, maximum(Δτ_sub_lin), length=200)

    # 1. Fit for Linear (Match the RoyalBlue of the dots)
    plot!(p1, Δτ_range, fit_lin.(Δτ_range);
        linestyle = :solid, 
        linewidth = 1.5,
        alpha = 0.4,       # Slightly transparent so it doesn't hide the points
        color = :royalblue, 
        label = "") 

    plot!(p1, Δτ_range, fit_lin_real.(Δτ_range);
        linestyle = :solid, 
        linewidth = 1.5,
        alpha = 0.4,       # Slightly transparent so it doesn't hide the points
        color = :midnightblue, 
        label = "")

    # 2. Fit for Quadratic (Match the Chocolate of the diamonds)
    plot!(p1, Δτ_range, fit_quad.(Δτ_range);
        linestyle = :solid, 
        linewidth = 2,
        alpha = 0.6,
        color = :chocolate, 
        label = "")
    plot!(p1, Δτ_range, fit_quad_real.(Δτ_range);
        linestyle = :solid,
        linewidth = 2,
        alpha = 0.6,
        color = :saddlebrown,
        label = "") 

    # Mark extrapolated values at Δτ = 0
    E0_lin  = fit_lin(0.0)
    E0_quad = fit_quad(0.0)

    E0_lin_real  = fit_lin_real(0.0)
    E0_quad_real = fit_quad_real(0.0)

    scatter!(p1, [0.0], [E0_lin];
         marker = :star5, markersize = 12, color = :royalblue,
         label = L"$E_0\mathrm{(lin)} = %$(round(E0_lin, digits=5))$")
    
    scatter!(p1, [0.0], [E0_lin_real];
         marker = :star5, markersize = 12, color = :midnightblue,
         label = L"$E_0\mathrm{(lin\ real)} = %$(round(E0_lin_real, digits=5))$")

    scatter!(p1, [0.0], [E0_quad];
         marker = :star5, markersize = 12, color = :chocolate,
         label = L"$E_0\mathrm{(quad)} = %$(round(E0_quad, digits=5))$")

    scatter!(p1, [0.0], [E0_quad_real];
            marker = :star5, markersize = 12, color = :saddlebrown,
            label = L"$E_0\mathrm{(quad\ real)} = %$(round(E0_quad_real, digits=5))$")

end
vmc_val = round(E_vmc / nr0_sq_32, digits=4)
vmc_err = round(Error_E_vmc / nr0_sq_32, digits=4)

# Fix: Ensure all three lines are normalized by nr0_sq_32
hline!(p1, [(E_vmc - Error_E_vmc)/nr0_sq_32, (E_vmc/nr0_sq_32), (E_vmc + Error_E_vmc)/nr0_sq_32];
       label = [L"$E_{\mathrm{VMC}} = %$vmc_val \pm %$vmc_err$" "" ""], 
       linestyle = :dot, 
       linewidth = 1.2,
       color = :black,
       alpha = 0.6)

savefig(p1, joinpath(@__DIR__, "..", "data", "plots", "DMC", "dmc_E_vs_dtau_N$(num_part)_nr0sq$(nr0_sq).pdf"))
savefig(p1, joinpath(@__DIR__, "..", "data", "plots", "DMC", "dmc_E_vs_dtau_N$(num_part)_nr0sq$(nr0_sq).png"))
println("Saved Plot 1")
display(p1)

# ------------------------------------------
# Plot 2: Energy vs num_walkers for each Δτ
# ------------------------------------------

dmc_path_walkers = joinpath(@__DIR__, "..", "data", "results", "DMC",
              "dmc_N$(num_part)_nr0sq$(nr0_sq)_quadratic.txt")
dmc_data_walkers = readdlm(dmc_path_walkers, '\t', Float64, skipstart=1)

num_walkers_vals_walkers = dmc_data_walkers[:, 1]
Δτ_vals_walkers = dmc_data_walkers[:, 2]
E_dmc_walkers = dmc_data_walkers[:, 3]
Error_E_dmc_walkers = dmc_data_walkers[:, 4]

unique_walkers = sort(unique(num_walkers_vals_walkers))
unique_Δτ      = sort(unique(Δτ_vals_walkers))
colors2 = palette(:plasma, length(unique_Δτ))

p2 = plot(
    xlabel  = L"1/N_{\mathrm{walkers}}",
    ylabel  = L"E/N \cdot (nr_0^2)^{-3/2}",
    title   = L"DMC\ \mathrm{Energy\ vs}\ N_{\mathrm{walkers}},\ N=30,\ nr_0^2=16",
    legend  = :topright,
    framestyle = :box,
    grid    = true,
    gridalpha = 0.3,
    xscale  = :log10,
    dpi = 600
)

unique_Δτ = [5*1e-5]

for (idx, dt) in enumerate(unique_Δτ)
    mask = Δτ_vals_walkers .== dt
    nw_sub  = num_walkers_vals_walkers[mask]
    E_sub   = E_dmc_walkers[mask]
    Err_sub = Error_E_dmc_walkers[mask]

    plot!(p2, 1 ./nw_sub, E_sub/nr0_sq_32,
          yerror   = Err_sub/nr0_sq_32,
          marker   = :circle,
          markersize = 5,
          linewidth  = 1.2,
          color    = colors2[idx])
end

savefig(p2, joinpath(@__DIR__, "..", "data", "plots", "DMC", "dmc_E_vs_walkers_N$(num_part)_nr0sq$(nr0_sq).pdf"))
savefig(p2, joinpath(@__DIR__, "..", "data", "plots", "DMC", "dmc_E_vs_walkers_N$(num_part)_nr0sq$(nr0_sq).png"))
println("Saved Plot 2")
display(p2)

