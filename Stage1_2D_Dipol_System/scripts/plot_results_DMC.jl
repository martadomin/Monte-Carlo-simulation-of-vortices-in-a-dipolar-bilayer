using DelimitedFiles, Plots, LaTeXStrings

ENV["PATH"] = "C:\\Users\\marta\\AppData\\Local\\Programs\\MiKTeX\\miktex\\bin\\x64;" * ENV["PATH"]
pgfplotsx()

# Load data
dmc_path = joinpath(@__DIR__, "..", "data", "results", "DMC",
           "dmc_N30_nr0sq16.0.txt")
dmc_data = readdlm(dmc_path, '\t', Float64, skipstart=1)

num_walkers_vals = dmc_data[:, 1]
Δτ_vals          = dmc_data[:, 2]
E_dmc            = dmc_data[:, 3]
Error_E_dmc      = dmc_data[:, 4]

num_part  = 30
nr0_sq    = 16.0
nr0_sq_32 = num_part * nr0_sq^(3/2)
E_tail    = 2π / sqrt(num_part)

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

savefig(p1, joinpath(@__DIR__, "..", "data", "plots", "DMC", "dmc_E_vs_dtau_N30_nr0sq16.pdf"))
println("Saved Plot 1")
display(p1)

# ------------------------------------------
# Plot 2: Energy vs num_walkers for each Δτ
# ------------------------------------------
colors2 = palette(:plasma, length(unique_Δτ))

p2 = plot(
    xlabel  = L"N_{\mathrm{walkers}}",
    ylabel  = L"E/N \cdot (nr_0^2)^{-3/2}",
    title   = L"DMC\ \mathrm{Energy\ vs}\ N_{\mathrm{walkers}},\ N=30,\ nr_0^2=16",
    legend  = :topright,
    framestyle = :box,
    grid    = true,
    gridalpha = 0.3
)

for (idx, dt) in enumerate(unique_Δτ)
    mask = Δτ_vals .== dt
    nw_sub  = num_walkers_vals[mask]
    E_sub   = E_norm[mask]
    Err_sub = Error_norm[mask]
    
    plot!(p2, nw_sub, E_sub,
          yerror   = Err_sub,
          label    = L"\Delta\tau = %$(dt)",
          marker   = :circle,
          markersize = 5,
          linewidth  = 2,
          color    = colors2[idx])
end

savefig(p2, joinpath(@__DIR__, "..", "data", "plots", "dmc_E_vs_walkers_N30_nr0sq16.pdf"))
println("Saved Plot 2")
display(p2)