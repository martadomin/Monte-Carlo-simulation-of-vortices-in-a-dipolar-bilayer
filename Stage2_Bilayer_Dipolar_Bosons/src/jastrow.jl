using Bessels
include(normpath(joinpath(@__DIR__, "shooting_method.jl")))

"""
    calculate_constants(L::Float64, R_match::Float64) -> Tuple{Float64, Float64, Float64}

Calculates the constants C1, C2, C3 for the two-body Jastrow factor of a 2D dipolar 
Bose gas. The Jastrow factor is defined piecewise: a short-range part proportional to 
K_0(2/√r) for r < R_match, and a long-range Reatto-Chester phononic form for r > R_match.

# Input:
- `L::Float64`: Length of the periodic simulation box.
- `R_match::Float64`: Matching distance separating the short- and long-range regimes. 
                      Must satisfy 0 < R_match < L/2.

# Output:
- `Tuple{Float64, Float64, Float64}`: The constants (C1, C2, C3) entering the Jastrow 
                                      factor. Must be computed in order C3 → C2 → C1.

# Notes:
The three constants are uniquely determined by:
  1. Continuity of f_2 at R_match
  2. Continuity of f_2' at R_match  
  3. Normalization f_2(L/2) = 1
See Ref. [Astrakharchik et al., PRL 2007] for details.
"""
function calculate_constants(L::Float64, R_match::Float64)
    C3 = besselk(1, 2/sqrt(R_match))/besselk(0, 2/sqrt(R_match)) * R_match^(-3/2) * (R_match^(-2) - (L-R_match)^(-2))^(-1)
    C2 = exp(4*C3/L)
    C1 = (C2 * exp(-C3/R_match) * exp(-C3/(L-R_match)))/besselk(0, 2/sqrt(R_match))
    return C1, C2, C3
end

# ================================================================
# SAME LAYER (AA and BB) — analytical Jastrow
# ================================================================

"""
    u_AA(r, R_match, L, Constants) -> Float64

Logarithm of the intra-layer two-body Jastrow factor, u_AA(r) = log(f_AA(r)).
Used in the Metropolis acceptance ratio and log-wavefunction evaluation.
Since f_BB = f_AA, this function is used for both AA and BB pairs.

Piecewise definition:
  - r < R_match : short-range, exact two-body scattering solution ~ K_0(2/√r)
  - r ≥ R_match : long-range, phononic Reatto-Chester form
  - r ≥ L/2    : cutoff, returns 0 (f_AA = 1)

# Input:
- `r::Float64`       : pair separation.
- `R_match::Float64` : matching radius between short and long range regimes.
- `L::Float64`       : simulation box length.
- `Constants`        : tuple (C1, C2, C3) from calculate_constants.
"""
function u_AA(r::Float64, R_match::Float64, L::Float64, 
              Constants::Tuple{Float64, Float64, Float64})::Float64
    if r < 1e-10    return 0.0 end
    if r >= L/2     return 0.0 end
    C1, C2, C3 = Constants
    if r < R_match
        return log(C1) + log(besselk(0, 2/sqrt(r)))
    else
        return log(C2) - C3/r - C3/(L-r)
    end
end

"""
    u_AA_prime(r, R_match, L, Constants) -> Float64

First derivative of the intra-layer log-Jastrow factor:
    u'_AA(r) = d/dr log(f_AA) = f'_AA / f_AA

Used in the drift force and kinetic energy calculation.
Since f_BB = f_AA, this function is used for both AA and BB pairs.

# Input:
- `r::Float64`       : pair separation.
- `R_match::Float64` : matching radius between short and long range regimes.
- `L::Float64`       : simulation box length.
- `Constants`        : tuple (C1, C2, C3) from calculate_constants.
"""
function u_AA_prime(r::Float64, R_match::Float64, L::Float64, 
                    Constants::Tuple{Float64, Float64, Float64})::Float64
    if r < 1e-10    return 0.0 end
    if r >= L/2     return 0.0 end
    _, _, C3 = Constants
    if r < R_match
        return besselk(1, 2/sqrt(r)) / besselk(0, 2/sqrt(r)) * r^(-3/2)
    else
        return C3/r^2 - C3/(L-r)^2
    end
end

"""
    u_AA_second(r, R_match, L, Constants) -> Float64

Second derivative of the intra-layer log-Jastrow factor:
    u''_AA(r) = d²/dr² log(f_AA) = f''_AA/f_AA - (f'_AA/f_AA)²

Used in the Laplacian term of the local kinetic energy.
Since f_BB = f_AA, this function is used for both AA and BB pairs.

# Input:
- `r::Float64`       : pair separation.
- `R_match::Float64` : matching radius between short and long range regimes.
- `L::Float64`       : simulation box length.
- `Constants`        : tuple (C1, C2, C3) from calculate_constants.
"""
function u_AA_second(r::Float64, R_match::Float64, L::Float64, 
                     Constants::Tuple{Float64, Float64, Float64})::Float64
    if r < 1e-10    return 0.0 end
    if r >= L/2     return 0.0 end
    _, _, C3 = Constants
    if r < R_match
        K0 = besselk(0, 2/sqrt(r))
        K1 = besselk(1, 2/sqrt(r))
        K2 = besselk(2, 2/sqrt(r))
        return (r)^(-3) * (((K0 * (K0 + K2))/2) - K1^2)/K0^2 - (3/2) * (K1/K0) * r^(-5/2)
    else
        return -2*C3/r^3 - 2*C3/(L-r)^3
    end
end

# ================================================================
# INTER LAYER (AB) — numerical Jastrow from shooting method
# ================================================================

"""
    u_AB(r, R0, r_grid, psi) -> Float64

Logarithm of the inter-layer two-body Jastrow factor, u_AB(r) = log(f_AB(r)).
Used in the Metropolis acceptance ratio and log-wavefunction evaluation.

f_AB is obtained numerically from the shooting method solution of the 
interlayer two-body Schrödinger equation. It is precomputed on a grid 
[r_min, R0] and evaluated here via linear interpolation.

Piecewise definition:
  - r ≤ r_grid[1] : below r_min, use first grid point
  - r_grid[1] < r < R0 : interpolated from precomputed grid
  - r ≥ R0        : f_AB = 1, returns log(1) = 0

# Input:
- `r::Float64`            : pair separation.
- `R0::Float64`           : cutoff radius, variational parameter.
- `r_grid::Vector{Float64}`: radial grid from shooting method.
- `psi::Vector{Float64}`  : f_AB values on r_grid from shooting method.
"""
function u_AB(r::Float64, R0::Float64, r_grid::Vector{Float64}, psi::Vector{Float64})::Float64
    # Ensure numerical stability for very small r
    if r < 1e-10 || r <= r_grid[1]
        return log(psi[1])
    end 

    if r >= R0 return 0.0 end
    
    i = searchsortedlast(r_grid, r)
    
    # Boundary check: if at the end of grid, return last value
    if i >= length(r_grid)
        return log(psi[end])
    end
    
    # Linear interpolation
    t = (r - r_grid[i]) / (r_grid[i+1] - r_grid[i])
    return (1-t)*log(psi[i]) + t*log(psi[i+1])
end

"""
    u_AB_prime(r, R0, r_grid, u_prime) -> Float64

First derivative of the inter-layer log-Jastrow factor:
    u'_AB(r) = d/dr log(f_AB) = f'_AB / f_AB

Evaluated via linear interpolation on the precomputed grid from 
the shooting method. Used in the drift force and kinetic energy.

# Input:
- `r::Float64`              : pair separation.
- `R0::Float64`             : cutoff radius, variational parameter.
- `r_grid::Vector{Float64}` : radial grid from shooting method.
- `u_prime::Vector{Float64}`: u'_AB values on r_grid from shooting method.
"""
function u_AB_prime(r::Float64, R0::Float64, r_grid::Vector{Float64}, 
                    u_prime::Vector{Float64})::Float64
    if r < 1e-10 || r <= r_grid[1]  return u_prime[1] end
    if r >= R0                       return 0.0 end
    i = searchsortedlast(r_grid, r)

        # Boundary check: if at the end of grid, return last value
    if i >= length(r_grid)
        return u_prime[end]
    end

    t = (r - r_grid[i]) / (r_grid[i+1] - r_grid[i])
    return (1-t)*u_prime[i] + t*u_prime[i+1]
end

"""
    u_AB_second(r, R0, r_grid, u_doubleprime) -> Float64

Second derivative of the inter-layer log-Jastrow factor:
    u''_AB(r) = d²/dr² log(f_AB) = f''_AB/f_AB - (f'_AB/f_AB)²

Evaluated via linear interpolation on the precomputed grid from 
the shooting method. Used in the Laplacian term of the local kinetic energy.

# Input:
- `r::Float64`                   : pair separation.
- `R0::Float64`                  : cutoff radius, variational parameter.
- `r_grid::Vector{Float64}`      : radial grid from shooting method.
- `u_doubleprime::Vector{Float64}`: u''_AB values on r_grid from shooting method.
"""
function u_AB_second(r::Float64, R0::Float64, r_grid::Vector{Float64}, 
                     u_doubleprime::Vector{Float64})::Float64
    if r < 1e-10 || r <= r_grid[1]  return u_doubleprime[1] end
    if r >= R0                       return 0.0 end
    i = searchsortedlast(r_grid, r)

    # Boundary check: if at the end of grid, return last value
    if i >= length(r_grid)
        return u_doubleprime[end]
    end

    t = (r - r_grid[i]) / (r_grid[i+1] - r_grid[i])
    return (1-t)*u_doubleprime[i] + t*u_doubleprime[i+1]
end