include("utils.jl")
include("jastrow.jl")

"""
    local_interaction_energy(x_coord, y_coord, L) -> Float64

Calculates the total pairwise dipole interaction energy for a set of particles in a
two-dimensional periodic box.

# Input:
- `x_coord::Vector{Float64}`: x-coordinates of all particles.
- `y_coord::Vector{Float64}`: y-coordinates of all particles.
- `L::Float64`: Length of the periodic box.

# Output:
- `Float64`: The total interaction energy, summed over all unique pairs.

# Notes
The interaction is given by `|r_i - r_j|^(-3)` for each pair `(i < j)`.
"""
function local_interaction_energy(x_A::Vector{Float64}, y_A::Vector{Float64}, x_B::Vector{Float64}, y_B::Vector{Float64}, L::Float64, h::Float64)::Float64
    E_int = 0.0
    N_half = length(x_A)
    h_sq = h^2

    @assert length(y_A) == N_half && length(x_B) == N_half && length(y_B) == N_half "Input coordinate vectors must have the same length"

    # AA Pairs =============
    @inbounds for i in 1:N_half
         @inbounds for j in (i + 1):N_half
             dx = get_periodic_difference(x_A[i], x_A[j], L)
             dy = get_periodic_difference(y_A[i], y_A[j], L)
             r_ij = dx^2 + dy^2
             if r_ij <= (L/2)^2
                 E_int += (r_ij)^(-3/2)
             end
         end
    end

    # BB Pairs =============
    @inbounds for α in 1:N_half
         @inbounds for β in (α + 1):N_half
             dx = get_periodic_difference(x_B[α], x_B[β], L)
             dy = get_periodic_difference(y_B[α], y_B[β], L)
             r_αβ = dx^2 + dy^2
             if r_αβ <= (L/2)^2
                 E_int += (r_αβ)^(-3/2)
             end
         end
    end

    # AB Pairs =============
    @inbounds for i in 1:N_half
        @inbounds for α in 1:N_half
            dx = get_periodic_difference(x_A[i], x_B[α], L)
            dy = get_periodic_difference(y_A[i], y_B[α], L)
            r2_iα = dx^2 + dy^2
            r3D_sq = r2_iα + h_sq
            if r2_iα <= (L/2)^2
                E_int += (r2_iα - 2*h_sq) / (r2_iα + h_sq)^(5/2)
            end
        end
    end

    return E_int
end

function phase_external_potential(x_A::Vector{Float64}, y_A::Vector{Float64},
                                 x_B::Vector{Float64}, y_B::Vector{Float64},
                                 x_vortex_A::Float64, y_vortex_A::Float64,
                                 x_vortex_B::Float64, y_vortex_B::Float64,
                                 L::Float64, lA::Float64, lB::Float64, rho2_min::Float64=1e-16)::Float64
    N_half = length(x_A)
    @assert length(y_A) == N_half && length(x_B) == N_half && length(y_B) == N_half "Input coordinate vectors must have the same length"

    V_ext = 0.0

    for i in 1:N_half

        dx_A = get_periodic_difference(x_A[i], x_vortex_A, L)
        dy_A = get_periodic_difference(y_A[i], y_vortex_A, L)
        dx_B = get_periodic_difference(x_B[i], x_vortex_B, L)
        dy_B = get_periodic_difference(y_B[i], y_vortex_B, L)

        rho2_A = max(dx_A^2 + dy_A^2, rho2_min)
        rho2_B = max(dx_B^2 + dy_B^2, rho2_min)

        V_ext += 0.5 * (lA^2) / (rho2_A) + 0.5 * (lB^2) / (rho2_B)        
    end
    return V_ext
end

function energy_estimators(x_A::Vector{Float64}, y_A::Vector{Float64},
                            x_B::Vector{Float64}, y_B::Vector{Float64},
                            x_vortex_A::Float64, y_vortex_A::Float64,
                            x_vortex_B::Float64, y_vortex_B::Float64,
                            L::Float64, h::Float64, lA::Float64, lB::Float64,
                            R_match::Float64,
                            Constants::Tuple{Float64, Float64, Float64},
                            R0::Float64,
                            itp_up,
                            itp_upp)::Tuple{Vector{Float64}, Vector{Float64}, Float64, Float64, Float64, Float64, Float64}

    N_half        = length(x_A)
    drift_x         = zeros(Float64, 2*N_half)
    drift_y         = zeros(Float64, 2*N_half)
    E_kin           = 0.0
    F_drift_sq      = 0.0
    Scalar_term_sum = 0.0
    E_int           = 0.0
    E_phase         = 0.0

    # ── Layer A particles ─────────────────────────────────────────
    @inbounds for k in 1:N_half
        F_x         = 0.0
        F_y         = 0.0
        scalar_term = 0.0

        @inbounds for j in 1:N_half

            # AA pairs
            if j != k
                dx = get_periodic_difference(x_A[k], x_A[j], L)
                dy = get_periodic_difference(y_A[k], y_A[j], L)
                r  = sqrt(dx^2 + dy^2)
                if r > 1e-10
                    du_dr   = u_AA_prime(r, R_match, L, Constants)
                    d2u_dr2 = u_AA_second(r, R_match, L, Constants)
                    F_x         += du_dr * (dx / r)
                    F_y         += du_dr * (dy / r)
                    scalar_term += d2u_dr2 + du_dr / r
                end
            end

            # AB pairs — same loop index j used as α
            dx = get_periodic_difference(x_A[k], x_B[j], L)
            dy = get_periodic_difference(y_A[k], y_B[j], L)
            r  = sqrt(dx^2 + dy^2)
            if r > 1e-10
                du_dr   = u_AB_prime(r, R0, itp_up)
                d2u_dr2 = u_AB_second(r, R0, itp_upp)
                F_x         += du_dr * (dx / r)
                F_y         += du_dr * (dy / r)
                scalar_term += d2u_dr2 + du_dr / r
            end

        end  # end j loop

        # Vortex contribution
        dx = get_periodic_difference(x_A[k], x_vortex_A, L)
        dy = get_periodic_difference(y_A[k], y_vortex_A, L)
        r  = sqrt(dx^2 + dy^2)
        if r > 1e-10
            du_dr   = u_vortex_prime(r, lA, L)
            d2u_dr2 = u_vortex_second(r, lA, L)
            F_x         += du_dr * (dx / r)
            F_y         += du_dr * (dy / r)
            scalar_term += d2u_dr2 + du_dr / r
        end

        drift_x[k]      = F_x
        drift_y[k]      = F_y
        F_drift_sq      += F_x^2 + F_y^2
        Scalar_term_sum += scalar_term
        E_kin           += F_x^2 + F_y^2 + scalar_term
    end

    # ── Layer B particles ─────────────────────────────────────────
    @inbounds for γ in 1:N_half
        F_x         = 0.0
        F_y         = 0.0
        scalar_term = 0.0

        @inbounds for β in 1:N_half

            # BB pairs
            if β != γ
                dx = get_periodic_difference(x_B[γ], x_B[β], L)
                dy = get_periodic_difference(y_B[γ], y_B[β], L)
                r  = sqrt(dx^2 + dy^2)
                if r > 1e-10
                    du_dr   = u_AA_prime(r, R_match, L, Constants)
                    d2u_dr2 = u_AA_second(r, R_match, L, Constants)
                    F_x         += du_dr * (dx / r)
                    F_y         += du_dr * (dy / r)
                    scalar_term += d2u_dr2 + du_dr / r
                end
            end

            # AB pairs — same loop index β used as i
            dx = get_periodic_difference(x_B[γ], x_A[β], L)
            dy = get_periodic_difference(y_B[γ], y_A[β], L)
            r  = sqrt(dx^2 + dy^2)
            if r > 1e-10
                du_dr   = u_AB_prime(r, R0, itp_up)
                d2u_dr2 = u_AB_second(r, R0, itp_upp)
                F_x         += du_dr * (dx / r)
                F_y         += du_dr * (dy / r)
                scalar_term += d2u_dr2 + du_dr / r
            end

        end  # end β loop

        # Vortex contribution
        dx = get_periodic_difference(x_B[γ], x_vortex_B, L)
        dy = get_periodic_difference(y_B[γ], y_vortex_B, L)
        r  = sqrt(dx^2 + dy^2)
        if r > 1e-10
            du_dr   = u_vortex_prime(r, lB, L)
            d2u_dr2 = u_vortex_second(r, lB, L)
            F_x         += du_dr * (dx / r)
            F_y         += du_dr * (dy / r)
            scalar_term += d2u_dr2 + du_dr / r
        end

        drift_x[N_half + γ] = F_x
        drift_y[N_half + γ] = F_y
        F_drift_sq          += F_x^2 + F_y^2
        Scalar_term_sum     += scalar_term
        E_kin               += F_x^2 + F_y^2 + scalar_term
    end

    # Interaction energy
    E_int = local_interaction_energy(x_A, y_A, x_B, y_B, L, h)

    # External potential energy
    E_phase = phase_external_potential(x_A, y_A, x_B, y_B,
                                        x_vortex_A, y_vortex_A, x_vortex_B, y_vortex_B,
                                        L, lA, lB)

    # Three kinetic estimators
    E_kin_std       = -0.5  * E_kin
    E_kin_drift     =  0.5  * F_drift_sq
    E_kin_laplacian = -0.25 * Scalar_term_sum
    E_total         = E_kin_std + E_int + E_phase
    
    return drift_x, drift_y,
           E_total,
           E_kin_drift + E_int + E_phase,
           E_kin_laplacian + E_int + E_phase,
           E_kin_std,
           E_int
end

function tail_energy(nr0sq::Float64, num_part::Int, h::Float64)::Float64
    ## AA and BB contribution:
    E_tail_AA_BB = (π * nr0sq^(3/2)) / sqrt(num_part)
    # ## AB contribution:
    E_tail_AB = (num_part * π)/(8*(num_part/(4*nr0sq) + h^2)^(3/2))

    return E_tail_AA_BB + E_tail_AB
end