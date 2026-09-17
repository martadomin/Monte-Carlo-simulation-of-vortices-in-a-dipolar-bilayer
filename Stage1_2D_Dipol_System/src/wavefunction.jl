# Stage1_2D_Dipol_System/src/wavefunction.jl

struct SingleLayerTrial <: TrialWavefunction
    R_match::Float64
    Constants::Tuple{Float64,Float64,Float64}

    function SingleLayerTrial(; R_match, Constants)
        new(R_match, Constants)
    end
end

"""
    logpsi(trial::SingleLayerTrial, coords, L) -> Float64

log|Ψ_T| for the single-layer trial wavefunction: log Ψ = Σ_{i<j} u_AA(r_ij).
Translated from the original compute_logΨ; u2 → u_AA (identical function,
renamed once Stage2 needed to distinguish it from u_AB — see common's
jastrow_common.jl).
"""
function logpsi(trial::SingleLayerTrial, coords::NamedTuple, L::Float64)::Float64
    x, y = coords.A.x, coords.A.y
    num_part = length(x)
    logΨ = 0.0
    @inbounds for i in 1:num_part, j in (i+1):num_part
        dx = get_periodic_difference(x[i], x[j], L)
        dy = get_periodic_difference(y[i], y[j], L)
        r = sqrt(dx^2 + dy^2)
        logΨ += u_AA(r, trial.R_match, L, trial.Constants)
    end
    return logΨ
end

"""
    Δlogpsi(trial::SingleLayerTrial, coords, species, id, x_new, y_new, L) -> Float64

Change in log|Ψ_T| from moving particle `id` to (x_new, y_new). Only valid
for a single-particle move. Translated from the original compute_ΔlogΨ.
"""
function Δlogpsi(trial::SingleLayerTrial, coords::NamedTuple, species::Symbol,
                  id::Int, x_new::Float64, y_new::Float64, L::Float64)::Float64
    x, y = coords.A.x, coords.A.y
    num_part = length(x)
    ΔlogΨ = 0.0
    @inbounds for i in 1:num_part
        if i != id
            r_old = sqrt(get_periodic_difference(x[i], x[id], L)^2 +
                         get_periodic_difference(y[i], y[id], L)^2)
            r_new = sqrt(get_periodic_difference(x[i], x_new, L)^2 +
                         get_periodic_difference(y[i], y_new, L)^2)
            ΔlogΨ += u_AA(r_new, trial.R_match, L, trial.Constants) -
                     u_AA(r_old, trial.R_match, L, trial.Constants)
        end
    end
    return ΔlogΨ
end

"""
    energy_estimators(trial::SingleLayerTrial, coords, L) -> NamedTuple

Translated from Stage1's original energy_estimators. u2_first_derivative/
u2_second_derivative → u_AA_prime/u_AA_second (identical functions, see
common's jastrow_common.jl).
"""
function energy_estimators(trial::SingleLayerTrial, coords::NamedTuple, L::Float64)
    x, y = coords.A.x, coords.A.y
    num_part = length(x)
    drift_x = zeros(Float64, num_part)
    drift_y = zeros(Float64, num_part)
    E_kin = F_drift_sq = Scalar_term_sum = E_int = 0.0

    @inbounds for k in 1:num_part
        F_x = F_y = scalar_term = 0.0
        for i in 1:num_part
            if i != k
                dx = get_periodic_difference(x[k], x[i], L)
                dy = get_periodic_difference(y[k], y[i], L)
                r  = sqrt(dx^2 + dy^2)
                if r > 1e-4
                    du_dr   = u_AA_prime(r, trial.R_match, L, trial.Constants)
                    d2u_dr2 = u_AA_second(r, trial.R_match, L, trial.Constants)
                    F_x += du_dr * (dx/r); F_y += du_dr * (dy/r)
                    scalar_term += d2u_dr2 + du_dr/r
                    if i < k && r <= L/2
                        E_int += (dx^2 + dy^2)^(-3/2)
                    end
                end
            end
        end
        drift_x[k], drift_y[k] = F_x, F_y
        F_drift_sq += F_x^2 + F_y^2
        Scalar_term_sum += scalar_term
        E_kin += F_x^2 + F_y^2 + scalar_term
    end

    E_kin_std       = -0.5  * E_kin
    E_kin_drift     =  0.5  * F_drift_sq
    E_kin_laplacian = -0.25 * Scalar_term_sum

    return (; drift = (A = (x = drift_x, y = drift_y),),
              E_total_local     = E_kin_std + E_int,
              E_total_drift     = E_kin_drift + E_int,
              E_total_laplacian = E_kin_laplacian + E_int,
              E_kinetic         = E_kin_std,
              E_interaction     = E_int)
end

mutable struct SingleLayerObservables
    n_xy::Matrix{Float64}
    gr::Vector{Float64}
    n_samples::Int
end

init_observables(::SingleLayerTrial, num_bins::Int) =
    SingleLayerObservables(zeros(num_bins, num_bins), zeros(num_bins), 0)

function record_observables!(::SingleLayerTrial, obs::SingleLayerObservables, coords::NamedTuple, L::Float64)
    accumulate_density!(obs.n_xy, coords.A.x, coords.A.y, L)
    accumulate_gr!(obs.gr, coords.A.x, coords.A.y, L)
    obs.n_samples += 1
end

"""
    tail_energy(num_part) -> Float64

Finite-size tail correction for the single-layer 2D dipolar Bose gas.
Single-layer analog of the bilayer tail_energy in
common/src/jastrow_interlayer.jl — no h dependence, no AB cross-term.
"""
function tail_energy(num_part::Int)::Float64
    return 2π / sqrt(num_part)
end