# common/src/plotting_diagnostics.jl
using Plots

"""
    plot_convergence(block_sizes, sigmas_std, sigmas_drift, sigmas_laplacian,
                      plateau_std, plateau_drift, plateau_laplacian; title="") -> Plots.Plot

Plots the blocking-error curves σ(B) vs. block size for all three energy
estimators (local/std, drift, laplacian) on one figure, each with its own
detected-plateau marker.
"""
function plot_convergence(block_sizes::Vector{Int},
                           sigmas_std::Vector{Float64}, sigmas_drift::Vector{Float64}, sigmas_laplacian::Vector{Float64},
                           plateau_std::Int, plateau_drift::Int, plateau_laplacian::Int;
                           title::String="")
    p = plot(block_sizes, sigmas_std; label="std", marker=:circle, color=:blue,
             xlabel="Block size", ylabel="σ", title=title, framestyle=:box)
    plot!(block_sizes, sigmas_drift; label="drift", marker=:square, color=:green)
    plot!(block_sizes, sigmas_laplacian; label="laplacian", marker=:diamond, color=:orange)
    vline!([plateau_std]; label="plateau (std)", linestyle=:dash, color=:blue)
    vline!([plateau_drift]; label="plateau (drift)", linestyle=:dash, color=:green)
    vline!([plateau_laplacian]; label="plateau (laplacian)", linestyle=:dash, color=:orange)
    return p
end

"""
    plot_dmc_trace(result, num_equil; title="") -> Plots.Plot

Plots the DMC energy trace E(step), with a vertical line at num_equil.
A genuine plateau should be visible before that line; if the trace is
still trending when it crosses num_equil, num_equil is too short.
"""
function plot_dmc_trace(result::NamedTuple, num_equil::Int; title::String="")
    p = plot(1:length(result.E_plot), result.E_plot; xlabel="DMC step", ylabel="E",
              title=title, label="E(step)", linewidth=1, framestyle=:box)
    vline!([num_equil]; label="num_equil", linestyle=:dash, color=:red)
    return p
end
