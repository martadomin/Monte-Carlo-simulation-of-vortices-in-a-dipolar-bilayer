using DelimitedFiles, Plots, LaTeXStrings, Polynomials, GLM, DataFrames

# At the top of the file, import explicitly:
import Polynomials: fit as poly_fit
import GLM: lm, coef, stderror

ENV["PATH"] = "C:\\Users\\marta\\AppData\\Local\\Programs\\MiKTeX\\miktex\\bin\\x64;" * ENV["PATH"]
gr()

default(
    fontfamily        = "Computer Modern",
    titlefontsize     = 11,
    guidefontsize     = 10,
    tickfontsize      = 8,
    legendfontsize    = 7,
    size              = (500, 380),
    linewidth         = 1.5,
    markersize        = 5,
    markerstrokewidth = 0.5,
    framestyle        = :box,
    grid              = true,
    gridalpha         = 0.25,
    gridlinewidth     = 0.5,
    gridstyle         = :dot,
    left_margin       = 5Plots.mm,
    bottom_margin     = 4Plots.mm,
    right_margin      = 3Plots.mm,
    top_margin        = 2Plots.mm,
    dpi               = 600
)

col_lin       = "#0072B2"
col_lin_real  = "#56B4E9"
col_quad      = "#D55E00"
col_quad_real = "#E69F00"

num_part    = 30
nr0_sq      = 32.0
nr0_sq_32   = num_part * nr0_sq^(3/2)
num_walkers = 200

function load_dmc(filename)
    path = joinpath(@__DIR__, "..", "data", "results", "DMC", filename)
    d = readdlm(path, '\t', Float64, skipstart=1)
    return d[:,1], d[:,2], d[:,3], d[:,4]
end

_, dt_lr, E_lr, err_lr = load_dmc("dmc_N$(num_part)_nr0sq$(nr0_sq)_linear.txt")
_, dt_qr, E_qr, err_qr = load_dmc("dmc_N$(num_part)_nr0sq$(nr0_sq)_quadratic.txt")

vmc_data  = readdlm(joinpath(@__DIR__, "..", "data", "results", "VMC",
                    "vmc_N$(num_part)_nr0sq$(nr0_sq).txt"), '\t', Float64, skipstart=1)
E_vmc     = vmc_data[1,5] / nr0_sq_32
err_vmc   = vmc_data[1,6] / nr0_sq_32

# Normalise
E_lr  ./= nr0_sq_32;  err_lr ./= nr0_sq_32
E_qr  ./= nr0_sq_32;  err_qr ./= nr0_sq_32

# Δτ mask
lo, hi = 1e-5, 1e-3
mask(dt) = (dt .>= lo) .& (dt .<= hi)
mlr = mask(dt_lr)
mqr = mask(dt_qr)

# # ── Helper: weighted GLM fit ──────────────────────────────────────────────────
# function glm_fit(dt, E, err, degree)
#     df = DataFrame(dt=dt, E=E, dt2=dt.^2, w=1.0./err.^2)
#     if degree == 1
#         m = lm(@formula(E ~ dt), df, wts=df.w)
#     else
#         m = lm(@formula(E ~ dt + dt2), df, wts=df.w)
#     end
#     return coef(m)[1], stderror(m)[1]   # E₀, σ_E₀
# end

# # Fits
# E0_l,  σ_l  = glm_fit(dt_l[ml],   E_l[ml],   err_l[ml],   1)
# E0_lr, σ_lr = glm_fit(dt_lr[mlr], E_lr[mlr], err_lr[mlr], 1)
# E0_q,  σ_q  = glm_fit(dt_q[mq],   E_q[mq],   err_q[mq],   2)
# E0_qr, σ_qr = glm_fit(dt_qr[mqr], E_qr[mqr], err_qr[mqr], 2)

# # Smooth fit curves from Polynomials.jl (for drawing only)
# fit_l_poly  = poly_fit(dt_l[ml],   E_l[ml],   1)
# fit_lr_poly = poly_fit(dt_lr[mlr], E_lr[mlr], 1)
# fit_q_poly  = poly_fit(dt_q[mq],   E_q[mq],   2)
# fit_qr_poly = poly_fit(dt_qr[mqr], E_qr[mqr], 2)

Δτ_range = range(0, hi, length=300)

# # ── Subplot 1: No branching ───────────────────────────────────────────────────
# p_nb = plot(
#     title  = L"\mathrm{No\ Branching}",
#     xlabel = L"\Delta\tau",
#     ylabel = L"E/N \cdot (nr_0^2)^{-3/2}",
#     legend = :bottomleft,
#     ylims  = (5.64, 5.675)
# )

# plot!(p_nb, dt_l[ml],  E_l[ml];  yerror=err_l[ml],  color=col_lin,  linestyle=:solid, marker=:circle,  label=L"\mathrm{Linear}")
# plot!(p_nb, dt_q[mq],  E_q[mq];  yerror=err_q[mq],  color=col_quad, linestyle=:solid, marker=:diamond, label=L"\mathrm{Quadratic}")
# plot!(p_nb, Δτ_range, fit_l_poly.(Δτ_range);  color=col_lin,  alpha=0.5, linewidth=1.5, label="")
# plot!(p_nb, Δτ_range, fit_q_poly.(Δτ_range);  color=col_quad, alpha=0.5, linewidth=1.5, label="")
# scatter!(p_nb, [0.0], [E0_l]; color=col_lin,  marker=:star5, markersize=12, label="")
# scatter!(p_nb, [0.0], [E0_q]; color=col_quad, marker=:star5, markersize=12, label="")
# hline!(p_nb, [E_vmc]; linestyle=:dot, linewidth=1.2, color=:black, alpha=0.7, label=L"E_\mathrm{VMC}")

# annotate!(p_nb, -1e-6, E0_l  - 0.004, text(L"E_0^{\rm lin}=%$(round(E0_l,  digits=5))",  :left, 6, col_lin))
# annotate!(p_nb, -1e-6, E0_q  + 0.003, text(L"E_0^{\rm quad}=%$(round(E0_q, digits=5))", :left, 6, col_quad))

# ── Subplot 2: With branching ─────────────────────────────────────────────────
p_b = plot(
    title  = L"\mathrm{With\ Branching}",
    xlabel = L"\Delta\tau",
    ylabel = "",
    legend = :bottomleft
    )

plot!(p_b, dt_lr[mlr], E_lr[mlr]; yerror=err_lr[mlr], color=col_lin_real,  linestyle=:solid, marker=:circle,  label=L"\mathrm{Linear}")
plot!(p_b, dt_qr[mqr], E_qr[mqr]; yerror=err_qr[mqr], color=col_quad_real, linestyle=:solid, marker=:diamond, label=L"\mathrm{Quadratic}")
# plot!(p_b, Δτ_range, fit_lr_poly.(Δτ_range); color=col_lin_real,  alpha=0.5, linewidth=1.5, label="")
# plot!(p_b, Δτ_range, fit_qr_poly.(Δτ_range); color=col_quad_real, alpha=0.5, linewidth=1.5, label="")
# scatter!(p_b, [0.0], [E0_lr]; color=col_lin_real,  marker=:star5, markersize=12, label="")
# scatter!(p_b, [0.0], [E0_qr]; color=col_quad_real, marker=:star5, markersize=12, label="")
# hline!(p_b, [E_vmc]; linestyle=:dot, linewidth=1.2, color=:black, alpha=0.7, label=L"E_\mathrm{VMC}")

# annotate!(p_b, 0, E0_lr - 0.015, text(L"E_0^{\rm lin}=%$(round(E0_lr,  digits=5))",  :left, 6, col_lin_real))
# annotate!(p_b, 0, E0_qr + 0.006, text(L"E_0^{\rm quad}=%$(round(E0_qr, digits=5))", :left, 6, col_quad_real))

# # ── Combine ───────────────────────────────────────────────────────────────────
# p_final = plot(p_nb, p_b,


#     layout      = (1, 2),
#     size        = (900, 400),
#     plot_title  = L"N=30,\quad nr_0^2=16,\quad N_w=200"
# )

savefig(p_b, joinpath(@__DIR__, "..", "data", "plots", "DMC",
        "dmc_E_vs_dtau_N$(num_part)_nr0sq$(nr0_sq).pdf"))
savefig(p_b, joinpath(@__DIR__, "..", "data", "plots", "DMC",
        "dmc_E_vs_dtau_N$(num_part)_nr0sq$(nr0_sq).png"))
display(p_b)