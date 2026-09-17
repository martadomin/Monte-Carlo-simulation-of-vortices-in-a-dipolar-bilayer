# Stage3_Vortex_Excitations/src/wavefunction.jl

struct VortexBilayerTrial <: TrialWavefunction
    R_match::Float64
    Constants::Tuple{Float64,Float64,Float64}
    R0::Float64
    itp_u
    itp_up
    itp_upp
    h::Float64
    lA::Float64
    lB::Float64
    x_vortex_A::Float64
    y_vortex_A::Float64
    x_vortex_B::Float64
    y_vortex_B::Float64

    function VortexBilayerTrial(; R_match, Constants, R0, itp_u, itp_up, itp_upp,
                                   h, lA, lB, x_vortex_A, y_vortex_A, x_vortex_B, y_vortex_B)
        new(R_match, Constants, R0, itp_u, itp_up, itp_upp,
            h, lA, lB, x_vortex_A, y_vortex_A, x_vortex_B, y_vortex_B)
    end
end

# ── Vortex one-body phase term — Stage3-only, no Bessels dependency ──

function u_vortex(r::Float64, l::Float64, L::Float64)::Float64
    r = max(r, 1e-8)
    if r > L/2    return 0.0 end
    return log(sin(pi * r / L)) * abs(l)
end

function u_vortex_prime(r::Float64, l::Float64, L::Float64)::Float64
    r = max(r, 1e-8)
    if r > L/2    return 0.0 end
    return (pi / L) * abs(l) * cot(pi * r / L)
end

function u_vortex_second(r::Float64, l::Float64, L::Float64)::Float64
    r = max(r, 1e-8)
    if r > L/2    return 0.0 end
    return - (pi^2 / L^2) * abs(l) * csc(pi * r / L)^2
end

# ── logpsi / Δlogpsi ──────────────────────────────────────────────
# Translated from the original compute_logΨ/compute_ΔlogΨ — Stage2's
# AA/BB/AB structure, plus the one-body vortex term per species.

function logpsi(trial::VortexBilayerTrial, coords::NamedTuple, L::Float64)::Float64
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

    @inbounds for i in 1:N_half
        rA = sqrt(get_periodic_difference(x_A[i], trial.x_vortex_A, L)^2 +
                  get_periodic_difference(y_A[i], trial.y_vortex_A, L)^2)
        logΨ += u_vortex(rA, trial.lA, L)
    end
    @inbounds for α in 1:N_half
        rB = sqrt(get_periodic_difference(x_B[α], trial.x_vortex_B, L)^2 +
                  get_periodic_difference(y_B[α], trial.y_vortex_B, L)^2)
        logΨ += u_vortex(rB, trial.lB, L)
    end
    return logΨ
end

function Δlogpsi(trial::VortexBilayerTrial, coords::NamedTuple, species::Symbol,
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
        rA_old = sqrt(get_periodic_difference(x_A[id], trial.x_vortex_A, L)^2 +
                      get_periodic_difference(y_A[id], trial.y_vortex_A, L)^2)
        rA_new = sqrt(get_periodic_difference(x_new, trial.x_vortex_A, L)^2 +
                      get_periodic_difference(y_new, trial.y_vortex_A, L)^2)
        ΔlogΨ += u_vortex(rA_new, trial.lA, L) - u_vortex(rA_old, trial.lA, L)
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
        rB_old = sqrt(get_periodic_difference(x_B[id], trial.x_vortex_B, L)^2 +
                      get_periodic_difference(y_B[id], trial.y_vortex_B, L)^2)
        rB_new = sqrt(get_periodic_difference(x_new, trial.x_vortex_B, L)^2 +
                      get_periodic_difference(y_new, trial.y_vortex_B, L)^2)
        ΔlogΨ += u_vortex(rB_new, trial.lB, L) - u_vortex(rB_old, trial.lB, L)
    end
    return ΔlogΨ
end

# ── energy_estimators ────────────────────────────────────────────

function phase_external_potential(x_A::Vector{Float64}, y_A::Vector{Float64},
                                   x_B::Vector{Float64}, y_B::Vector{Float64},
                                   x_vortex_A::Float64, y_vortex_A::Float64,
                                   x_vortex_B::Float64, y_vortex_B::Float64,
                                   L::Float64, lA::Float64, lB::Float64,
                                   rho2_min::Float64=1e-16)::Float64
    N_half = length(x_A)
    V_ext = 0.0
    @inbounds for i in 1:N_half
        rho2_A = max(get_periodic_difference(x_A[i], x_vortex_A, L)^2 +
                     get_periodic_difference(y_A[i], y_vortex_A, L)^2, rho2_min)
        rho2_B = max(get_periodic_difference(x_B[i], x_vortex_B, L)^2 +
                     get_periodic_difference(y_B[i], y_vortex_B, L)^2, rho2_min)
        V_ext += 0.5 * lA^2 / rho2_A + 0.5 * lB^2 / rho2_B
    end
    return V_ext
end

function energy_estimators(trial::VortexBilayerTrial, coords::NamedTuple, L::Float64)
    x_A, y_A = coords.A.x, coords.A.y
    x_B, y_B = coords.B.x, coords.B.y
    N_half   = length(x_A)

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
                    F_x += du_dr * (dx/r); F_y += du_dr * (dy/r)
                    scalar_term += d2u_dr2 + du_dr/r
                end
            end
            dx = get_periodic_difference(x_A[k], x_B[j], L)
            dy = get_periodic_difference(y_A[k], y_B[j], L)
            r  = sqrt(dx^2 + dy^2)
            if r > 1e-10
                du_dr   = u_AB_prime(r, trial.R0, trial.itp_up)
                d2u_dr2 = u_AB_second(r, trial.R0, trial.itp_upp)
                F_x += du_dr * (dx/r); F_y += du_dr * (dy/r)
                scalar_term += d2u_dr2 + du_dr/r
            end
        end
        dx = get_periodic_difference(x_A[k], trial.x_vortex_A, L)
        dy = get_periodic_difference(y_A[k], trial.y_vortex_A, L)
        r  = sqrt(dx^2 + dy^2)
        if r > 1e-10
            du_dr   = u_vortex_prime(r, trial.lA, L)
            d2u_dr2 = u_vortex_second(r, trial.lA, L)
            F_x += du_dr * (dx/r); F_y += du_dr * (dy/r)
            scalar_term += d2u_dr2 + du_dr/r
        end
        drift_x[k], drift_y[k] = F_x, F_y
        F_drift_sq += F_x^2 + F_y^2
        Scalar_term_sum += scalar_term
        E_kin += F_x^2 + F_y^2 + scalar_term
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
                    F_x += du_dr * (dx/r); F_y += du_dr * (dy/r)
                    scalar_term += d2u_dr2 + du_dr/r
                end
            end
            dx = get_periodic_difference(x_B[γ], x_A[β], L)
            dy = get_periodic_difference(y_B[γ], y_A[β], L)
            r  = sqrt(dx^2 + dy^2)
            if r > 1e-10
                du_dr   = u_AB_prime(r, trial.R0, trial.itp_up)
                d2u_dr2 = u_AB_second(r, trial.R0, trial.itp_upp)
                F_x += du_dr * (dx/r); F_y += du_dr * (dy/r)
                scalar_term += d2u_dr2 + du_dr/r
            end
        end
        dx = get_periodic_difference(x_B[γ], trial.x_vortex_B, L)
        dy = get_periodic_difference(y_B[γ], trial.y_vortex_B, L)
        r  = sqrt(dx^2 + dy^2)
        if r > 1e-10
            du_dr   = u_vortex_prime(r, trial.lB, L)
            d2u_dr2 = u_vortex_second(r, trial.lB, L)
            F_x += du_dr * (dx/r); F_y += du_dr * (dy/r)
            scalar_term += d2u_dr2 + du_dr/r
        end
        drift_x[N_half+γ], drift_y[N_half+γ] = F_x, F_y
        F_drift_sq += F_x^2 + F_y^2
        Scalar_term_sum += scalar_term
        E_kin += F_x^2 + F_y^2 + scalar_term
    end

    E_interaction = local_interaction_energy(x_A, y_A, x_B, y_B, L, trial.h)
    E_phase = phase_external_potential(x_A, y_A, x_B, y_B,
                  trial.x_vortex_A, trial.y_vortex_A, trial.x_vortex_B, trial.y_vortex_B,
                  L, trial.lA, trial.lB)

    E_kin_std       = -0.5  * E_kin
    E_kin_drift     =  0.5  * F_drift_sq
    E_kin_laplacian = -0.25 * Scalar_term_sum

    return (; drift = (A = (x = drift_x[1:N_half], y = drift_y[1:N_half]),
                        B = (x = drift_x[N_half+1:end], y = drift_y[N_half+1:end])),
              E_total_local     = E_kin_std + E_interaction + E_phase,
              E_total_drift     = E_kin_drift + E_interaction + E_phase,
              E_total_laplacian = E_kin_laplacian + E_interaction + E_phase,
              E_kinetic         = E_kin_std,
              E_interaction     = E_interaction)
end

# ── observables — own copy, structurally identical to Stage2's ────

mutable struct VortexBilayerObservables
    n_xy_A::Matrix{Float64}
    n_xy_B::Matrix{Float64}
    gAA_r::Vector{Float64}
    gBB_r::Vector{Float64}
    gAB_r::Vector{Float64}
    n_samples::Int
end

init_observables(::VortexBilayerTrial, num_bins::Int) = VortexBilayerObservables(
    zeros(num_bins, num_bins), zeros(num_bins, num_bins),
    zeros(num_bins), zeros(num_bins), zeros(num_bins), 0)

function record_observables!(::VortexBilayerTrial, obs::VortexBilayerObservables, coords::NamedTuple, L::Float64)
    accumulate_density!(obs.n_xy_A, coords.A.x, coords.A.y, L)
    accumulate_density!(obs.n_xy_B, coords.B.x, coords.B.y, L)
    accumulate_gr!(obs.gAA_r, coords.A.x, coords.A.y, L)
    accumulate_gr!(obs.gBB_r, coords.B.x, coords.B.y, L)
    accumulate_g_AB_r!(obs.gAB_r, coords.A.x, coords.A.y, coords.B.x, coords.B.y, L)
    obs.n_samples += 1
end