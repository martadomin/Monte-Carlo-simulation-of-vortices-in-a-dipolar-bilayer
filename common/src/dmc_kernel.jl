# common/src/dmc_kernel.jl

using Random

"""
    diffusion_step(coords, L, D, Δτ) -> NamedTuple

Proposes new walker positions via free diffusion (no drift), with periodic
wrapping into the box, for every species in `coords`. One step of the
diffusion half of a DMC move.
"""
function diffusion_step(coords::NamedTuple, L::Float64, D::Float64, Δτ::Float64)::NamedTuple
    return map(coords) do c
        (x = wrap_position.(c.x .+ sqrt(2*D*Δτ) .* randn(length(c.x)), L),
         y = wrap_position.(c.y .+ sqrt(2*D*Δτ) .* randn(length(c.y)), L))
    end
end

"""
    weight_update(E_loc_old, E_loc_new, E_ref, Δτ) -> Float64

Importance-sampling weight multiplier for one DMC step, using the
trapezoidal-rule approximation of the branching Green's function.
Returns 0.0 (walker dies) if either local energy is NaN. The exponent
is clamped to [-50, log(10)] to avoid overflow/underflow in `exp`.
"""
function weight_update(E_loc_old::Float64, E_loc_new::Float64,
                       E_ref::Float64, Δτ::Float64)::Float64
    if isnan(E_loc_old) || isnan(E_loc_new)
        return 0.0
    end
    exponent = -Δτ * (0.5*(E_loc_old + E_loc_new) - E_ref)
    return exp(clamp(exponent, -50.0, log(10.0)))
end

"""
    branching_step(weights) -> Vector{Int}

Converts walker weights into integer copy counts via stochastic rounding
(floor(w + rand())). NaN weights are treated as 0 (walker dies); weights
are capped at 3 to bound the maximum number of copies from a single walker
in one step.
"""
function branching_step(weights::Vector{Float64})::Vector{Int}
    safe_weights = [isnan(w) ? 0.0 : clamp(w, 0.0, 3.0) for w in weights]
    return [floor(Int, w + rand()) for w in safe_weights]
end

"""
    population_control(num_walkers, num_target, avg_E_loc, Δτ) -> Float64

Updates the reference energy E_ref to push the walker population back
toward num_target, using proportional feedback with gain α = 1.0.
"""
function population_control(num_walkers::Int, num_target::Int,
                             avg_E_loc::Float64, Δτ::Float64)::Float64
    α = 1.0
    return avg_E_loc - (α / Δτ) * log(num_walkers / num_target)
end