using DelimitedFiles, Plots, LaTeXStrings, Polynomials
import Polynomials: fit as poly_fit

ENV["PATH"] = "C:\\Users\\marta\\AppData\\Local\\Programs\\MiKTeX\\miktex\\bin\\x64;" * ENV["PATH"]
gr()

default(
    fontfamily="Computer Modern", titlefontsize=11, guidefontsize=10,
    tickfontsize=8, legendfontsize=7, size=(500,380), linewidth=1.5,
    markersize=5, markerstrokewidth=0.5, framestyle=:box,
    grid=true, gridalpha=0.25, gridlinewidth=0.5, gridstyle=:dot,
    left_margin=5Plots.mm, bottom_margin=4Plots.mm,
    right_margin=3Plots.mm, top_margin=2Plots.mm, dpi=600
)

col_quad_real = "#E69F00"

num_part    = 30
nr0_sq      = 1024.0
nr0_sq_32   = num_part * nr0_sq^(3/2)
num_walkers = 100

ref_data    = readdlm(joinpath(@__DIR__, "..", "data", "sweep_results",
              "delta_optimal_N$(num_part)_DMC.txt"), '\t', Float64, skipstart=1)
idx_ref     = findfirst(x -> isapprox(x, nr0_sq; atol=1e-8), ref_data[:,1])
Δτ_ref      = idx_ref === nothing ? 0.0 : ref_data[idx_ref, 4]

function load_dmc(filename)
    d = readdlm(joinpath(@__DIR__, "..", "data", "results", "DMC", filename),
                Float64, skipstart=1)
    return d[:,1], d[:,2], d[:,3], d[:,4]
end

_, dt_qr, E_qr, err_qr = load_dmc(
    "dmc_N$(num_part)_nr0sq$(nr0_sq)_Nw$(num_walkers)_quadratic.txt")

vmc_data = readdlm(joinpath(@__DIR__, "..", "data", "results", "VMC",
           "vmc_N$(num_part)_nr0sq$(nr0_sq).txt"), '\t', Float64, skipstart=1)
E_vmc    = vmc_data[1,5] / nr0_sq_32

E_qr   ./= nr0_sq_32
err_qr ./= nr0_sq_32

mqr = (dt_qr .>= 1e-7) .& (dt_qr .<= 1e-3)

# poly_fit for both curve and star — guaranteed consistent
fit_qr_poly = poly_fit(dt_qr[mqr], E_qr[mqr], 2)
E0_qr       = fit_qr_poly(0.0)

Δτ_range = range(0, maximum(dt_qr[mqr]), length=300)

p_b = plot(
    title  = L"\mathrm{With\ Branching}",
    xlabel = L"\Delta\tau",
    ylabel = L"E/N \cdot (nr_0^2)^{-3/2}"
)

plot!(p_b, dt_qr[mqr], E_qr[mqr];
    yerror=err_qr[mqr], color=col_quad_real,
    linestyle=:solid, marker=:diamond, label=L"\mathrm{Quadratic}")

plot!(p_b, Δτ_range, fit_qr_poly.(Δτ_range);
    color=col_quad_real, alpha=0.5, linewidth=1.5, label="")

vline!(p_b, [Δτ_ref];
    color=:black, linestyle=:dash,
    label=L"\Delta\tau_{\mathrm{ref}}", linewidth=1.0)

scatter!(p_b, [0.0], [E0_qr];
    color=col_quad_real, marker=:star5, markersize=12,
    label=L"E_0 = %$(round(E0_qr, digits=5))")

savefig(p_b, joinpath(@__DIR__, "..", "data", "plots", "DMC",
        "dmc_E_vs_dtau_N$(num_part)_nr0sq$(nr0_sq).pdf"))
savefig(p_b, joinpath(@__DIR__, "..", "data", "plots", "DMC",
        "dmc_E_vs_dtau_N$(num_part)_nr0sq$(nr0_sq).png"))
display(p_b)