# common/src/jastrow_interlayer.jl
#
# Inter-layer Jastrow factor u_AB, for two-species (bilayer) systems only
# — Stage2 and Stage3, not Stage1. f_AB is obtained numerically via the
# shooting method (see shooting_method.jl) and consumed here purely as an
# opaque interpolant object, so this file has no direct dependency on
# shooting_method.jl at load time.

"""
    u_AB(r, R0, itp_u) -> Float64

Logarithm of the inter-layer two-body Jastrow factor. r ≥ R0 → f_AB = 1,
i.e. returns 0.
"""
function u_AB(r::Float64, R0::Float64, itp_u)::Float64
    r >= R0 && return 0.0
    return itp_u(r)
end

"""
    u_AB_prime(r, R0, itp_up) -> Float64

First derivative of the inter-layer log-Jastrow factor.
"""
function u_AB_prime(r::Float64, R0::Float64, itp_up)::Float64
   r >= R0 && return 0.0
    return itp_up(r)
end

"""
    u_AB_second(r, R0, itp_upp) -> Float64

Second derivative of the inter-layer log-Jastrow factor.
"""
function u_AB_second(r::Float64, R0::Float64, itp_upp)::Float64
    r >= R0 && return 0.0
    return itp_upp(r)
end

"""
    local_interaction_energy(x_A, y_A, x_B, y_B, L, h) -> Float64

Total pairwise dipole interaction energy (AA + BB + AB) for a bilayer
system. Identical between Stage2 (bilayer) and Stage3 (vortex bilayer) —
the interaction Hamiltonian itself doesn't change when vortices are added,
only the trial wavefunction and the external potential do.
"""
function local_interaction_energy(x_A::Vector{Float64}, y_A::Vector{Float64},
                                   x_B::Vector{Float64}, y_B::Vector{Float64},
                                   L::Float64, h::Float64)::Float64
    E_int = 0.0
    N_half = length(x_A)
    h_sq = h^2
    @assert length(y_A) == N_half && length(x_B) == N_half && length(y_B) == N_half "Input coordinate vectors must have the same length"

    @inbounds for i in 1:N_half, j in (i+1):N_half
        dx = get_periodic_difference(x_A[i], x_A[j], L)
        dy = get_periodic_difference(y_A[i], y_A[j], L)
        r_ij = dx^2 + dy^2
        r_ij <= (L/2)^2 && (E_int += r_ij^(-3/2))
    end
    @inbounds for α in 1:N_half, β in (α+1):N_half
        dx = get_periodic_difference(x_B[α], x_B[β], L)
        dy = get_periodic_difference(y_B[α], y_B[β], L)
        r_αβ = dx^2 + dy^2
        r_αβ <= (L/2)^2 && (E_int += r_αβ^(-3/2))
    end
    @inbounds for i in 1:N_half, α in 1:N_half
        dx = get_periodic_difference(x_A[i], x_B[α], L)
        dy = get_periodic_difference(y_A[i], y_B[α], L)
        r2_iα = dx^2 + dy^2
        r2_iα <= (L/2)^2 && (E_int += (r2_iα - 2*h_sq) / (r2_iα + h_sq)^(5/2))
    end
    return E_int
end

"""
    tail_energy(nr0sq, num_part, h) -> Float64

Finite-size tail correction for the bilayer dipolar interaction energy,
accounting for pair correlations beyond the periodic box's cutoff at L/2.
A post-hoc correction applied to a blocked energy estimate — not called
by energy_estimators itself. Applies to Stage2/Stage3; not used by
Stage1's single-layer system.
"""
function tail_energy(nr0sq::Float64, num_part::Int, h::Float64)::Float64
    E_tail_AA_BB = (π * nr0sq^(3/2)) / sqrt(num_part)
    E_tail_AB    = (num_part * π) / (8*(num_part/(4*nr0sq) + h^2)^(3/2))
    return E_tail_AA_BB + E_tail_AB
end