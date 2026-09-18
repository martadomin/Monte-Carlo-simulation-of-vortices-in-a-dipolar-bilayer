# common/src/data_management.jl

using JLD2, Dates, Plots

"""
    result_path(stage_dir, kind, params) -> String

Builds a deterministic path under `<stage_dir>/data/results/<kind>/`,
encoding `params` in the filename. One function, used by every stage's
scripts, so the naming convention and folder layout live in one place
instead of being hand-built separately in each script.
"""
function result_path(stage_dir::String, kind::String, params::NamedTuple)::String
    param_str = join(["$(k)$(v)" for (k,v) in pairs(params)], "_")
    return joinpath(stage_dir, "data", "results", kind, "$(kind)_$(param_str).jld2")
end

"""
    save_run(path, params, result; overwrite=false)

Saves one run's raw result (a NamedTuple, e.g. from metropolis()/dmc())
to `path`, alongside its parameters and a timestamp.

- `overwrite=true`  → path is replaced entirely, starting fresh.
- `overwrite=false` (default) → if `path` already has runs saved, this one
  is added as a new entry rather than replacing what's there — useful for
  rerunning the same parameters to accumulate statistics or check
  reproducibility.
"""
function save_run(path::String, params::NamedTuple, result::NamedTuple; overwrite::Bool=false)
    mkpath(dirname(path))
    entry = (; params, result, timestamp=now())

    if overwrite || !isfile(path)
        jldsave(path; run_1=entry)
        return
    end

    n = jldopen(f -> length(keys(f)), path, "r")
    jldopen(path, "r+") do f
        f["run_$(n+1)"] = entry
    end
end

"""
    load_run(path; run=nothing)

Loads run(s) saved by save_run. With `run` omitted, returns every run in
the file as a Vector, in order. With `run` given (an Int), returns just
that one run's (; params, result, timestamp).
"""
function load_run(path::String; run::Union{Int,Nothing}=nothing)
    jldopen(path, "r") do f
        ks = sort(collect(keys(f)); by = k -> parse(Int, replace(k, "run_" => "")))
        return run === nothing ? [f[k] for k in ks] : f["run_$run"]
    end
end

"""
    plot_convergence(block_sizes, sigmas_std, sigmas_drift, sigmas_laplacian,
                      plateau_std, plateau_drift, plateau_laplacian; title="") -> Plots.Plot

Plots the blocking-error curves σ(B) vs. block size for all three energy
estimators (local/std, drift, laplacian) on one figure, each with its own
detected-plateau marker. The three curves should broadly agree on where
they plateau — a visible mismatch between them is worth investigating,
since all three are estimating the same physical quantity.
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