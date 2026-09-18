# Stage2_Bilayer_Dipolar_Bosons/src/h_sweep_analysis.jl
#
# Shared data assembly for the h-dependence figures (Fig1, Fig1_inset).
# Kept separate from optimization.jl since this is post-hoc analysis of
# already-saved results, not something that runs a simulation.

using Interpolations

"""
    load_h_sweep(N, nr0sq, h_vals, stage_dir) -> (h_found, E_corr_vals, err_vals)

Loads each h's saved R0_optimum result, adds the tail correction, and
returns only the h values that were actually found (skipping any missing
with a warning) — so callers never need to handle NaN/missing themselves.
"""
function load_h_sweep(N::Int, nr0sq::Float64, h_vals::Vector{Float64}, stage_dir::String)
    h_found, E_corr_vals, err_vals = Float64[], Float64[], Float64[]
    for h in h_vals
        path = result_path(stage_dir, "R0_optimum", (N=N, nr0sq=nr0sq, h=h))
        if !isfile(path)
            @warn "No saved R0_optimum for h=$h — skipping"
            continue
        end
        r = load_run(path; run=1).result
        push!(h_found, h)
        push!(E_corr_vals, r.E_opt + tail_energy(nr0sq, N, h))
        push!(err_vals, r.err_opt)
    end
    return h_found, E_corr_vals, err_vals
end

"""
    single_layer_reference(N, nr0sq, stage1_dir) -> Float64

Loads Stage1's saved optimum for a single layer at density nr0sq/2 with
N÷2 particles — the h→∞ decoupled-layers limit — and returns its
tail-corrected E/N. Returns NaN (with a warning) if not found.
"""
function single_layer_reference(N::Int, nr0sq::Float64, stage1_dir::String)::Float64
    N_half, nr0sq_half = N ÷ 2, nr0sq / 2
    path = result_path(stage1_dir, "Rmatch_optimum", (N=N_half, L=sqrt(N_half / nr0sq_half)))
    if !isfile(path)
        @warn "No saved single-layer reference at $path"
        return NaN
    end
    r = load_run(path; run=1).result
    return r.E_opt_total / N_half + tail_energy(N_half)
end

"""
    exact_binding_energy_interp(stage_dir) -> Function

Loads the exact dimer binding energy table and returns a linear
interpolant ε_b(h) — call it directly rather than hand-rolling the
interpolation each time.
"""
function exact_binding_energy_interp(stage_dir::String)
    data = readdlm(joinpath(stage_dir, "data", "binding_energy_dimer", "dimer_binding_energy.txt"),
                    '\t', skipstart=1)
    h_exact, eb_exact = Float64.(data[:, 1]), Float64.(data[:, 2])
    return LinearInterpolation(h_exact, eb_exact; extrapolation_bc=Flat())
end

"""
    load_h_sweep_dmc(N, nr0sq, h_vals, stage_dir) -> (h_found, E_corr_vals, err_vals)

Same as load_h_sweep, but reads the extrapolated DMC result (kind="DMC",
saved by num_walkers_convergence.jl) instead of the VMC optimum — the
tail correction is added here too, exactly as for VMC, since it wasn't
part of the DMC extrapolation itself.
"""
function load_h_sweep_dmc(N::Int, nr0sq::Float64, h_vals::Vector{Float64}, stage_dir::String)
    h_found, E_corr_vals, err_vals = Float64[], Float64[], Float64[]
    for h in h_vals
        path = result_path(stage_dir, "DMC", (N=N, nr0sq=nr0sq, h=h))
        if !isfile(path)
            @warn "No saved DMC result for h=$h — skipping"
            continue
        end
        r = load_run(path; run=1).result
        push!(h_found, h)
        push!(E_corr_vals, r.E_extrapolated + tail_energy(nr0sq, N, h))
        push!(err_vals, r.E_extrapolated_err)
    end
    return h_found, E_corr_vals, err_vals
end