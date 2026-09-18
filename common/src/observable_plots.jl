# common/src/observable_plots.jl

using Plots

"""
    plot_density_map(n_xy, xy_bins; title="") -> Plots.Plot

Heatmap of a normalized 2D density n(x,y). Call normalize_density! on the
raw accumulated n_xy first — this function only plots.
"""
function plot_density_map(n_xy::Matrix{Float64}, xy_bins::Vector{Float64}; title::String="")
    heatmap(xy_bins, xy_bins, n_xy'; xlabel="x", ylabel="y", title=title,
            color=:viridis, aspect_ratio=:equal, framestyle=:box)
end

"""
    plot_gr(r_vals, g_r; label="", title="") -> Plots.Plot
"""
function plot_gr(r_vals::Vector{Float64}, g_r::Vector{Float64}; label::String="", title::String="")
    plot(r_vals, g_r; xlabel="r", ylabel="g(r)", label=label, title=title,
         linewidth=2, framestyle=:box)
end

"""
    extrapolate_mixed(O_VMC, O_DMC; guard=1e-12) -> (extr_linear, extr_quadratic)

DMC extrapolated estimators for operators that don't commute with H (so
the mixed estimator alone is biased, unlike for energy):
    extr_linear    = 2·O_DMC - O_VMC
    extr_quadratic = O_DMC² / O_VMC   (guaranteed non-negative)
`guard` keeps the ratio form well-defined where O_VMC is zero or
near-zero (e.g. deep inside a vortex core, or short-range g(r) bins
suppressed by the Jastrow factor).
"""
function extrapolate_mixed(O_VMC, O_DMC; guard=1e-12)
    size(O_VMC) == size(O_DMC) || error("The VMC and DMC observables must have the same size.")
    extr_linear = 2 .* O_DMC .- O_VMC
    O_vmc_safe = sign.(O_VMC) .* max.(abs.(O_VMC), guard)
    extr_quadratic = O_DMC .^ 2 ./ O_vmc_safe
    return extr_linear, extr_quadratic
end