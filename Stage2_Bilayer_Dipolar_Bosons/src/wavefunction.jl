# Stage2_Bilayer_Dipolar_Bosons/src/wavefunction.jl

struct BilayerTrial <: TrialWavefunction
    R_match::Float64
    Constants::Tuple{Float64,Float64,Float64}
    R0::Float64
    itp_u
    itp_up
    itp_upp
    h::Float64

    function BilayerTrial(; R_match, Constants, R0, itp_u, itp_up, itp_upp, h)
        new(R_match, Constants, R0, itp_u, itp_up, itp_upp, h)
    end
end

"""
    logpsi(trial::BilayerTrial, coords, L) -> Float64

log|Ψ_T| = Σ_{i<j} u_AA(r_ij) [layer A] + Σ_{α<β} u_AA(r_αβ) [layer B]
         + Σ_{i,α} u_AB(r_iα)
Translated from the original compute_logΨ.
"""
function logpsi(trial::BilayerTrial, coords::NamedTuple, L::Float64)::Float64
    x_A, y_A = coords.A.x, coords.A.y
    x_B, y_B = coords.B.x, coords.B.y
    N_half = length(x_A)
    logΨ = 0.0

    @inbounds for i in 1:N_half, j in (i+1):N_half
        dx = get_periodic_difference(x_A[i], x_A[j], L)
        dy = get_periodic_difference(y_A[i], y_A[j], L)
        logΨ += u_AA(sqrt(dx^2+dy^2), trial.R_match, L, trial.Constants)
    end
    @inbounds for α in 1:N_half, β in (α+1):N_half
        dx = get_periodic_difference(x_B[α], x_B[β], L)
        dy = get_periodic_difference(y_B[α], y_B[β], L)
        logΨ += u_AA(sqrt(dx^2+dy^2), trial.R_match, L, trial.Constants)
    end
    @inbounds for i in 1:N_half, α in 1:N_half
        dx = get_periodic_difference(x_A[i], x_B[α], L)
        dy = get_periodic_difference(y_A[i], y_B[α], L)
        logΨ += u_AB(sqrt(dx^2+dy^2), trial.R0, trial.itp_u)
    end
    return logΨ
end

"""
    Δlogpsi(trial::BilayerTrial, coords, species, id, x_new, y_new, L) -> Float64

Translated from the original compute_ΔlogΨ.
"""
function Δlogpsi(trial::BilayerTrial, coords::NamedTuple, species::Symbol,
                  id::Int, x_new::Float64, y_new::Float64, L::Float64)::Float64
    x_A, y_A = coords.A.x, coords.A.y
    x_B, y_B = coords.B.x, coords.B.y
    N_half = length(x_A)
    ΔlogΨ = 0.0

    if species == :A
        @inbounds for j in 1:N_half
            if j != id
                r_old = sqrt(get_periodic_difference(x_A[id], x_A[j], L)^2 +
                             get_periodic_difference(y_A[id], y_A[j], L)^2)
                r_new = sqrt(get_periodic_difference(x_new, x_A[j], L)^2 +
                             get_periodic_difference(y_new, y_A[j], L)^2)
                ΔlogΨ += u_AA(r_new, trial.R_match, L, trial.Constants) -
                         u_AA(r_old, trial.R_match, L, trial.Constants)
            end
            r_old = sqrt(get_periodic_difference(x_A[id], x_B[j], L)^2 +
                         get_periodic_difference(y_A[id], y_B[j], L)^2)
            r_new = sqrt(get_periodic_difference(x_new, x_B[j], L)^2 +
                         get_periodic_difference(y_new, y_B[j], L)^2)
            ΔlogΨ += u_AB(r_new, trial.R0, trial.itp_u) - u_AB(r_old, trial.R0, trial.itp_u)
        end
    else  # :B
        @inbounds for j in 1:N_half
            if j != id
                r_old = sqrt(get_periodic_difference(x_B[id], x_B[j], L)^2 +
                             get_periodic_difference(y_B[id], y_B[j], L)^2)
                r_new = sqrt(get_periodic_difference(x_new, x_B[j], L)^2 +
                             get_periodic_difference(y_new, y_B[j], L)^2)
                ΔlogΨ += u_AA(r_new, trial.R_match, L, trial.Constants) -
                         u_AA(r_old, trial.R_match, L, trial.Constants)
            end
            r_old = sqrt(get_periodic_difference(x_B[id], x_A[j], L)^2 +
                         get_periodic_difference(y_B[id], y_A[j], L)^2)
            r_new = sqrt(get_periodic_difference(x_new, x_A[j], L)^2 +
                         get_periodic_difference(y_new, y_A[j], L)^2)
            ΔlogΨ += u_AB(r_new, trial.R0, trial.itp_u) - u_AB(r_old, trial.R0, trial.itp_u)
        end
    end
    return ΔlogΨ
end

"""
    energy_estimators(trial::BilayerTrial, coords, L) -> NamedTuple

Translated from the original energy_estimators. Calls common's
local_interaction_energy rather than duplicating the pairwise sum.
"""
function energy_estimators(trial::BilayerTrial, coords::NamedTuple, L::Float64)
    x_A, y_A = coords.A.x, coords.A.y
    x_B, y_B = coords.B.x, coords.B.y
    N_half = length(x_A)
    drift_x = zeros(Float64, 2*N_half)
    drift_y = zeros(Float64, 2*N_half)
    E_kin = F_drift_sq = Scalar_term_sum = 0.0

    @inbounds for k in 1:N_half
        F_x = F_y = scalar_term = 0.0
        @inbounds for j in 1:N_half
            if j != k
                dx = get_periodic_difference(x_A[k], x_A[j], L)
                dy = get_periodic_difference(y_A[k], y_A[j], L)
                r  = sqrt(dx^2 + dy^2)
                if r > 1e-10
                    du_dr   = u_AA_prime(r, trial.R_match, L, trial.Constants)
                    d2u_dr2 = u_AA_second(r, trial.R_match, L, trial.Constants)
                    F_x += du_dr*(dx/r); F_y += du_dr*(dy/r)
                    scalar_term += d2u_dr2 + du_dr/r
                end
            end
            dx = get_periodic_difference(x_A[k], x_B[j], L)
            dy = get_periodic_difference(y_A[k], y_B[j], L)
            r  = sqrt(dx^2 + dy^2)
            if r > 1e-10
                du_dr   = u_AB_prime(r, trial.R0, trial.itp_up)
                d2u_dr2 = u_AB_second(r, trial.R0, trial.itp_upp)
                F_x += du_dr*(dx/r); F_y += du_dr*(dy/r)
                scalar_term += d2u_dr2 + du_dr/r
            end
        end
        drift_x[k], drift_y[k] = F_x, F_y
        F_drift_sq += F_x^2+F_y^2; Scalar_term_sum += scalar_term
        E_kin += F_x^2+F_y^2+scalar_term
    end

    @inbounds for γ in 1:N_half
        F_x = F_y = scalar_term = 0.0
        @inbounds for β in 1:N_half
            if β != γ
                dx = get_periodic_difference(x_B[γ], x_B[β], L)
                dy = get_periodic_difference(y_B[γ], y_B[β], L)
                r  = sqrt(dx^2 + dy^2)
                if r > 1e-10
                    du_dr   = u_AA_prime(r, trial.R_match, L, trial.Constants)
                    d2u_dr2 = u_AA_second(r, trial.R_match, L, trial.Constants)
                    F_x += du_dr*(dx/r); F_y += du_dr*(dy/r)
                    scalar_term += d2u_dr2 + du_dr/r
                end
            end
            dx = get_periodic_difference(x_B[γ], x_A[β], L)
            dy = get_periodic_difference(y_B[γ], y_A[β], L)
            r  = sqrt(dx^2 + dy^2)
            if r > 1e-10
                du_dr   = u_AB_prime(r, trial.R0, trial.itp_up)
                d2u_dr2 = u_AB_second(r, trial.R0, trial.itp_upp)
                F_x += du_dr*(dx/r); F_y += du_dr*(dy/r)
                scalar_term += d2u_dr2 + du_dr/r
            end
        end
        drift_x[N_half+γ], drift_y[N_half+γ] = F_x, F_y
        F_drift_sq += F_x^2+F_y^2; Scalar_term_sum += scalar_term
        E_kin += F_x^2+F_y^2+scalar_term
    end

    E_interaction = local_interaction_energy(x_A, y_A, x_B, y_B, L, trial.h)
    E_kin_std       = -0.5  * E_kin
    E_kin_drift     =  0.5  * F_drift_sq
    E_kin_laplacian = -0.25 * Scalar_term_sum

    return (; drift = (A=(x=drift_x[1:N_half], y=drift_y[1:N_half]),
                        B=(x=drift_x[N_half+1:end], y=drift_y[N_half+1:end])),
              E_total_local     = E_kin_std + E_interaction,
              E_total_drift     = E_kin_drift + E_interaction,
              E_total_laplacian = E_kin_laplacian + E_interaction,
              E_kinetic         = E_kin_std,
              E_interaction     = E_interaction)
end

mutable struct BilayerObservables
    n_xy_A::Matrix{Float64}
    n_xy_B::Matrix{Float64}
    gAA_r::Vector{Float64}
    gBB_r::Vector{Float64}
    gAB_r::Vector{Float64}
    n_samples::Int
end

init_observables(::BilayerTrial, num_bins::Int) = BilayerObservables(
    zeros(num_bins, num_bins), zeros(num_bins, num_bins),
    zeros(num_bins), zeros(num_bins), zeros(num_bins), 0)

function record_observables!(::BilayerTrial, obs::BilayerObservables, coords::NamedTuple, L::Float64)
    accumulate_density!(obs.n_xy_A, coords.A.x, coords.A.y, L)
    accumulate_density!(obs.n_xy_B, coords.B.x, coords.B.y, L)
    accumulate_gr!(obs.gAA_r, coords.A.x, coords.A.y, L)
    accumulate_gr!(obs.gBB_r, coords.B.x, coords.B.y, L)
    accumulate_g_AB_r!(obs.gAB_r, coords.A.x, coords.A.y, coords.B.x, coords.B.y, L)
    obs.n_samples += 1
end