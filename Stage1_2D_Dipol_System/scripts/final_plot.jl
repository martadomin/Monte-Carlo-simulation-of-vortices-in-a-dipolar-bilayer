# ENV["PATH"] = "C:\\Users\\marta\\AppData\\Local\\Programs\\MiKTeX\\miktex\\bin\\x64;" * ENV["PATH"]
using DelimitedFiles, Plots, LaTeXStrings, PGFPlotsX

E_vmc = Float64[]
error_E_vmc = Float64[]

E_dmc = Float64[]
error_E_dmc = Float64[]

nr0_sq_vals = [16.0, 32.0, 48.0, 64.0, 96.0, 128.0, 196.0, 256.0]
num_part = 30
L = sqrt.(num_part ./ nr0_sq_vals)

for nr0_sq in nr0_sq_vals
        nr0_sq_32 = num_part * nr0_sq^(3/2)

        vmc_path = joinpath(@__DIR__, "..", "data", "results", "VMC",
                "vmc_N$(num_part)_nr0sq$(nr0_sq).txt")
        vmc_data = readdlm(vmc_path, '\t', Float64, skipstart=1)
        push!(E_vmc, vmc_data[1, 5] / nr0_sq_32)     
        push!(error_E_vmc, vmc_data[1, 6] / nr0_sq_32)

        dmc_path = joinpath(@__DIR__, "..", "data", "results", "DMC",
           "dmc_N$(num_part)_nr0sq$(nr0_sq)_quadratic_def.txt")
        dmc_data = readdlm(dmc_path, '\t', Float64, skipstart=1)
        push!(E_dmc, dmc_data[1, 4] / nr0_sq_32)
        push!(error_E_dmc, dmc_data[1, 5] / nr0_sq_32)
end

# Paper fit: E/N = a1*(nr0^2)^(3/2) + a2*(nr0^2)^(5/4) + a3*(nr0^2)^(1/2)
a1, a2, a3 = 4.536, 4.38, 1.2
nr0_range  = collect(LinRange(10.0, 300.0, 500))
E_paper    = @. a1*nr0_range^(3/2) + a2*nr0_range^(5/4) + a3*nr0_range^(1/2)
E_paper_normalized = E_paper ./ nr0_range.^(3/2)
E_tail = 2 * π / sqrt(num_part)

E_fit_data = @. a1*nr0_sq_vals^(3/2) + a2*nr0_sq_vals^(5/4) + a3*nr0_sq_vals^(1/2) - E_tail * nr0_sq_vals^(3/2)
E_fit_data_normalized = E_fit_data ./ nr0_sq_vals.^(3/2)

println(E_fit_data_normalized)

println("Relative error DMC vs paper fit: ", round.(abs.(E_dmc .- E_fit_data_normalized) ./ abs.(E_fit_data_normalized) * 100, sigdigits=3), "%")

# pgfplotsx()
gr()

xticks_vals = [0.0, 50.0, 100.0, 150.0, 200.0, 250.0, 300.0]
xticks_labels = [L"0", L"50", L"100", L"150", L"200", L"250", L"300"]
yticks_vals = [5.5, 6.0, 6.5, 7.0, 7.5]
yticks_labels = [L"5.5", L"6.0", L"6.5", L"7.0", L"7.5"]

p2 = plot(nr0_range, E_paper_normalized,
          label=L"\mathrm{DMC\ fit\ (Astrakharchik\ 2007)}",
          xlabel=L"nr_0^2",
          ylabel=L"E/N \cdot (nr_0^2)^{-3/2}",
          title=L"\mathrm{VMC\ vs\ DMC},\ N = %$(num_part)",
          linewidth=1,
          color=:black,
          framestyle=:box,
          legend=:topright,
          xticks=(xticks_vals, xticks_labels),
          yticks=(yticks_vals, yticks_labels))

# scatter!(nr0_sq_vals, E_vmc .+ E_tail,
#          yerror=error_E_vmc,
#          label=L"\mathrm{VMC,\ }N = %$(num_part)",
#          marker=:circle,
#          markersize=6,
#          color=:red)
println(E_dmc .+ E_tail)
println(E_vmc .+ E_tail)

nr0 = 16.0
println((a1*nr0^(3/2) + a2*nr0^(5/4) + a3*nr0^(1/2)) / nr0^(3/2))

scatter!(nr0_sq_vals, E_dmc .+ E_tail,
         yerror=error_E_dmc,
         label=L"\mathrm{DMC\ results,\ }N = %$(num_part)",
         marker=:square,
         markersize=6,
         color=:green)

savefig(joinpath(@__DIR__, "..", "data", "plots",
        "plot_vmc_vs_paper_N$(num_part).pdf"))
println("Saved final plot")
display(p2)

# # g(r) plot for different nr0^2 values
# nr0_sq_available = [96.0]
# colors = palette(:viridis, length(nr0_sq_available))

# gr()

# p_gr = plot(
#     xlabel = L"r/r_0",
#     ylabel = L"g_2(r)",
#     title  = L"g_2(r),\ N = %$(num_part)",
#     legend = :topright,
#     framestyle = :box,
#     grid = true,
#     gridalpha = 0.3,
#     size = (800, 500)
# )

# for (idx, nr0_sq) in enumerate(nr0_sq_available)
#     gr_path = joinpath(@__DIR__, "..", "data", "results",
#               "gr_N$(num_part)_nr0sq$(nr0_sq).txt")
#     if !isfile(gr_path)
#         println("Skipping g(r) for nr0^2 = $nr0_sq — file not found")
#         continue
#     end
#     gr_data = readdlm(gr_path, '\t', Float64, skipstart=1)
#     r_vals  = gr_data[:, 1]
#     gr      = gr_data[:, 2]
#     plot!(r_vals, gr,
#           label  = L"nr_0^2 = %$(Int(nr0_sq))",
#           color  = colors[idx],
#           linewidth = 2)
# end

# savefig(joinpath(@__DIR__, "..", "data", "plots",
#         "plot_gr_N$(num_part).pdf"))
# println("Saved g(r) plot")
# display(p_gr)