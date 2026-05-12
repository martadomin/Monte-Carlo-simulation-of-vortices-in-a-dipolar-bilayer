using DelimitedFiles, Plots, LaTeXStrings, Polynomials

ENV["PATH"] = "C:\\Users\\marta\\AppData\\Local\\Programs\\MiKTeX\\miktex\\bin\\x64;" * ENV["PATH"]
# pgfplotsx()
gr()

default(
    fontfamily        = "Computer Modern",
    titlefontsize     = 12,
    guidefontsize     = 11,
    tickfontsize      = 9,
    legendfontsize    = 7,
    annotationfontsize = 9,
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
    right_margin      = 5Plots.mm,
    top_margin        = 1Plots.mm,
    dpi               = 600
)

# ── Colorblind-safe palette (Wong 2011) ───────────────────────────────────────
col_lin      = "#0072B2"   # blue        → linear no branching
col_lin_real = "#56B4E9"   # sky blue    → linear real
col_quad     = "#D55E00"   # vermillion  → quadratic no branching
col_quad_real= "#E69F00"   # orange      → quadratic real

num_part  = 30
nr0_sq    = 16.0
nr0_sq_32 = num_part * nr0_sq^(3/2)
E_tail    = 2π / sqrt(num_part)
println("E_tail = ", E_tail)
num_walkers = 200

# ── Load data ─────────────────────────────────────────────────────────────────
function load_dmc(filename, num_part, nr0_sq)
    path = joinpath(@__DIR__, "..", "data", "results", "DMC", filename)
    d = readdlm(path, '\t', Float64, skipstart=1)
    return d[:, 1], d[:, 2], d[:, 3], d[:, 4]   # nw, Δτ, E, err
end

nw_lin,       dt_lin,       E_lin,       err_lin       = load_dmc("dmc_N$(num_part)_nr0sq$(nr0_sq)_Nw$(num_walkers)_linear_no_branching.txt",   num_part, nr0_sq)
nw_lin_real,  dt_lin_real,  E_lin_real,  err_lin_real  = load_dmc("dmc_N$(num_part)_nr0sq$(nr0_sq)_Nw$(num_walkers)_linear.txt",                num_part, nr0_sq)
nw_quad,      dt_quad,      E_quad,      err_quad      = load_dmc("dmc_N$(num_part)_nr0sq$(nr0_sq)_Nw$(num_walkers)_quadratic_no_branching.txt", num_part, nr0_sq)
nw_quad_real, dt_quad_real, E_quad_real, err_quad_real = load_dmc("dmc_N$(num_part)_nr0sq$(nr0_sq)_Nw$(num_walkers)_quadratic.txt",             num_part, nr0_sq)

vmc_data    = readdlm(joinpath(@__DIR__, "..", "data", "results", "VMC", "vmc_N$(num_part)_nr0sq$(nr0_sq).txt"), '\t', Float64, skipstart=1)
E_vmc       = vmc_data[1, 5] / nr0_sq_32
Error_E_vmc = vmc_data[1, 6] / nr0_sq_32

# ── Normalise ─────────────────────────────────────────────────────────────────
E_n_lin       = E_lin       ./ nr0_sq_32;  err_n_lin       = err_lin       ./ nr0_sq_32
E_n_lin_real  = E_lin_real  ./ nr0_sq_32;  err_n_lin_real  = err_lin_real  ./ nr0_sq_32
E_n_quad      = E_quad      ./ nr0_sq_32;  err_n_quad      = err_quad      ./ nr0_sq_32
E_n_quad_real = E_quad_real ./ nr0_sq_32;  err_n_quad_real = err_quad_real ./ nr0_sq_32

unique_walkers = sort(unique(nw_lin))

# ── Plot 1: E vs Δτ ───────────────────────────────────────────────────────────
p1 = plot(
    xlabel  = L"\Delta\tau",
    ylabel  = L"E/N \cdot (nr_0^2)^{-3/2}",
    title   = L"\mathrm{DMC\ Energy\ vs}\ \Delta\tau,\quad N=30,\quad nr_0^2=16",
    legend  = :bottomleft,
    legendcolumns = 1,
    ylims = (5.53, 5.68)
)

for nw in unique_walkers
    dt_range = (1e-5, 1e-3)
    function submask(nw_arr, dt_arr)
        (nw_arr .== nw) .& (dt_arr .>= dt_range[1]) .& (dt_arr .<= dt_range[2])
    end

    m_lin       = submask(nw_lin,       dt_lin)
    m_lin_real  = submask(nw_lin_real,  dt_lin_real)
    m_quad      = submask(nw_quad,      dt_quad)
    m_quad_real = submask(nw_quad_real, dt_quad_real)

    # ── Subsets ────────────────────────────────────────────────────────────────
    Δτ_l,  E_l,  err_l  = dt_lin[m_lin],             E_n_lin[m_lin],       err_n_lin[m_lin]
    Δτ_lr, E_lr, err_lr = dt_lin_real[m_lin_real],   E_n_lin_real[m_lin_real],  err_n_lin_real[m_lin_real]
    Δτ_q,  E_q,  err_q  = dt_quad[m_quad],           E_n_quad[m_quad],     err_n_quad[m_quad]
    Δτ_qr, E_qr, err_qr = dt_quad_real[m_quad_real], E_n_quad_real[m_quad_real], err_n_quad_real[m_quad_real]

    # ── Legend labels: only series type, NOT E₀ values ────────────────────────
    # (E₀ values go as annotations near Δτ = 0 below)
    lbl_l  = L"N_w=%$nw\ \mathrm{(lin,\ no\ branch)}"
    lbl_lr = L"N_w=%$nw\ \mathrm{(lin,\ with\ branch)}"
    lbl_q  = L"N_w=%$nw\ \mathrm{(quad,\ no\ branch)}"
    lbl_qr = L"N_w=%$nw\ \mathrm{(quad,\ with\ branch)}"

    # ── Data series ───────────────────────────────────────────────────────────
    plot!(p1, Δτ_l,  E_l;  yerror=err_l,  label=lbl_l,  color=col_lin,       linestyle=:dash,   marker=:circle,  markersize=5)
    plot!(p1, Δτ_lr, E_lr; yerror=err_lr, label=lbl_lr, color=col_lin_real,  linestyle=:solid,  marker=:circle,  markersize=5)
    plot!(p1, Δτ_q,  E_q;  yerror=err_q,  label=lbl_q,  color=col_quad,      linestyle=:dash,   marker=:diamond, markersize=6)
    plot!(p1, Δτ_qr, E_qr; yerror=err_qr, label=lbl_qr, color=col_quad_real, linestyle=:solid,  marker=:diamond, markersize=6)

    # ── Polynomial fits ───────────────────────────────────────────────────────
    Δτ_fit = range(0, maximum(Δτ_l), length=200)
    fit_l  = fit(Δτ_l,  E_l,  1)
    fit_lr = fit(Δτ_lr, E_lr, 1)
    fit_q  = fit(Δτ_q,  E_q,  2)
    fit_qr = fit(Δτ_qr, E_qr, 2)

    plot!(p1, Δτ_fit, fit_l.(Δτ_fit);  color=col_lin,       alpha=0.5, linewidth=1.5, linestyle=:solid, label="")
    plot!(p1, Δτ_fit, fit_lr.(Δτ_fit); color=col_lin_real,  alpha=0.5, linewidth=1.5, linestyle=:solid, label="")
    plot!(p1, Δτ_fit, fit_q.(Δτ_fit);  color=col_quad,      alpha=0.6, linewidth=2.0, linestyle=:solid, label="")
    plot!(p1, Δτ_fit, fit_qr.(Δτ_fit); color=col_quad_real, alpha=0.6, linewidth=2.0, linestyle=:solid, label="")

    # ── Extrapolated E₀ — scatter with NO label ───────────────────────────────
    E0_l  = fit_l(0.0);   E0_lr = fit_lr(0.0)
    E0_q  = fit_q(0.0);   E0_qr = fit_qr(0.0)

    scatter!(p1, [0.0], [E0_l];  color=col_lin,       marker=:star5, markersize=10, label="")
    scatter!(p1, [0.0], [E0_lr]; color=col_lin_real,  marker=:star5, markersize=10, label="")
    scatter!(p1, [0.0], [E0_q];  color=col_quad,      marker=:star5, markersize=10, label="")
    scatter!(p1, [0.0], [E0_qr]; color=col_quad_real, marker=:star5, markersize=10, label="")

    # ── Annotations: E₀ values placed right of the star markers ───────────────
 
    annotate!(p1, 0.00, 5.655,  text(L"E_0^\mathrm{l}=%$(round(E0_l,  digits=5))", :left, 7, col_lin))
    annotate!(p1, 0.00, 5.59, text(L"E_0^\mathrm{lr}=%$(round(E0_lr, digits=5))", :left, 7, col_lin_real))
    annotate!(p1, 0.00, 5.675,  text(L"E_0^\mathrm{q}=%$(round(E0_q,  digits=5))", :left, 7, col_quad))
    annotate!(p1, 0.00, 5.615, text(L"E_0^\mathrm{qr}=%$(round(E0_qr, digits=5))", :left, 7, col_quad_real))
end


# ── VMC reference: hline with NO legend entry, annotated at the right edge ───
hline!(p1, [E_vmc];
       linestyle=:dot, linewidth=1.5, color=:black, alpha=0.7, label="")

# Annotation placed at right edge, just above the line
x_right = maximum(dt_lin) * 0.98
annotate!(p1, x_right, E_vmc + (E_vmc * 0.0015),
    text(L"E_\mathrm{VMC}=%$(round(E_vmc,digits=4))\pm%$(round(Error_E_vmc,digits=4))",
         :right, 7, :black))

savefig(p1, joinpath(@__DIR__, "..", "data", "plots", "DMC",
        "dmc_E_vs_dtau_N$(num_part)_nr0sq$(nr0_sq).png")
        )
savefig(p1, joinpath(@__DIR__, "..", "data", "plots", "DMC",
        "dmc_E_vs_dtau_N$(num_part)_nr0sq$(nr0_sq).pdf"))
println("Saved Plot 1")
display(p1)

# ── Plot 2: E vs N_walkers ────────────────────────────────────────────────────
p2 = plot(
    xlabel  = L"1/N_{\mathrm{walkers}}",
    ylabel  = L"E/N \cdot (nr_0^2)^{-3/2}",
    title   = L"\mathrm{DMC\ Energy\ vs}\ N_{\mathrm{walkers}},\quad N=30,\quad nr_0^2=16",
    legend  = :topright,
    xscale  = :log10
)

for dt in [1e-4]
    mask   = dt_quad .== dt
    nw_sub = nw_quad[mask]
    E_sub  = E_n_quad[mask]
    err_sub= err_n_quad[mask]
    scatter!(p2, 1 ./ nw_sub, E_sub;
             yerror=err_sub, label=L"\Delta\tau=%$(dt)",
             marker=:circle, markersize=5, color=col_quad)
end

savefig(p2, joinpath(@__DIR__, "..", "data", "plots", "DMC",
        "dmc_E_vs_walkers_N$(num_part)_nr0sq$(nr0_sq).pdf"))
println("Saved Plot 2")
display(p2)