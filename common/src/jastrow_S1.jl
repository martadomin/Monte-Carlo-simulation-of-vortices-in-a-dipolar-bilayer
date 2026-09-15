using Bessels

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

"""
    u2(r::Float64) -> Float64
Defines the logarithm of the two-body wave function
"""
function u2(r::Float64, R_match::Float64, L::Float64, Constants::Tuple{Float64, Float64, Float64})::Float64
    # Ensure numerical stability for very small r
    if r < 1e-10
        return 0.0
    end

    if r>=L/2 return 0.0 end

    C1, C2, C3 = Constants
    if r < R_match
        return log(C1) + log(besselk(0, 2/sqrt(r)))
    else
        return log(C2) - C3/r - C3/(L-r)
    end
end

"""
    u2_first_derivative(r::Float64, R_match::Float64, L::Float64, Constants::Tuple{Float64, Float64, Float64})::Float64
Defines the first derivative of the logarithm of the two-body wave function
"""
function u2_first_derivative(r::Float64, R_match::Float64, L::Float64, Constants::Tuple{Float64, Float64, Float64})::Float64
    # Ensure numerical stability for very small r
    if r < 1e-10
        return 0.0
    end

    if r>=L/2 return 0.0 end

    _, _, C3 = Constants
    if r < R_match
        return besselk(1, 2/sqrt(r)) / besselk(0, 2/sqrt(r)) * r^(-3/2)
    else
        return C3/r^2 - C3/(L-r)^2
    end
end

"""
    u2_second_derivative(r::Float64, R_match::Float64, L::Float64, Constants::Tuple{Float64, Float64, Float64})::Float64
Defines the second derivative of the logarithm of the two-body wave function
"""
function u2_second_derivative(r::Float64, R_match::Float64, L::Float64, Constants::Tuple{Float64, Float64, Float64})::Float64
    # Ensure numerical stability for very small r
    if r < 1e-10
        return 0.0
    end

    if r>=L/2 return 0.0 end

    _, _, C3 = Constants
    if r < R_match
        K0 = besselk(0, 2/sqrt(r))
        K1 = besselk(1, 2/sqrt(r))
        K2 = besselk(2, 2/sqrt(r))
        return (r)^(-3) * (((K0 * (K0 + K2))/2) - K1^2)/K0^2 - (3/2) * (K1/K0) * r^(-5/2)
    else
        return - 2* C3/r^3 - 2 * C3/(L-r)^3
    end
end

