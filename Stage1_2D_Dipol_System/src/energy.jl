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
function local_interaction_energy(x_coord::Vector{Float64}, y_coord::Vector{Float64}, L::Float64)::Float64
    E_int = 0.0
    num_part = length(x_coord)

    @inbounds for i in 1:num_part
        @inbounds for j in (i + 1):num_part
            dx = get_periodic_difference(x_coord[i], x_coord[j], L)
            dy = get_periodic_difference(y_coord[i], y_coord[j], L)
            r = sqrt(dx^2 + dy^2)
            if r > L/2
                E_int += (dx^2 + dy^2)^(-3/2)
            else
                E_int += (dx^2 + dy^2)^(-3/2)
            end
        end
    end
    return E_int
end

"""
    local_kinetic_energy(xcoord, ycoord, L, R_match, Constants) -> Float64

Evaluates the local kinetic energy of the 2D dipolar Bose gas using the
logarithmic-derivative form of the Jastrow wavefunction.

# Input:
- `xcoord::Vector{Float64}`: x-coordinates of all particles.
- `ycoord::Vector{Float64}`: y-coordinates of all particles.
- `L::Float64`: Box length (periodic boundary conditions assumed).
- `R_match::Float64`: Matching distance between short- and long-range regimes of the Jastrow factor.
- `Constants::Tuple{Float64, Float64, Float64}`: Tuple (C1, C2, C3) of Jastrow constants.

# Output:
- `Float64`: Local kinetic energy in dimensionless units.
"""
function local_kinetic_energy(xcoord::Vector{Float64}, ycoord::Vector{Float64}, L::Float64, R_match::Float64, 
                                Constants::Tuple{Float64, Float64, Float64})::Tuple{Float64, Float64, Float64}
    
    num_part = length(xcoord)
    E_kin = 0.0
    F_drift_sq = 0.0
    Scalar_term_sum = 0.0

    @inbounds for k in 1:num_part
        F_drift_x = 0.0
        F_drift_y = 0.0
        scalar_term = 0.0
        for i in 1:num_part
            if i != k
                dx = get_periodic_difference(xcoord[k], xcoord[i], L)
                dy = get_periodic_difference(ycoord[k], ycoord[i], L)
                r = sqrt(dx^2 + dy^2)
                du_dr = u2_first_derivative(r, R_match, L, Constants)
                d2u_dr2 = u2_second_derivative(r, R_match, L, Constants)
                F_drift_x += du_dr * (dx / r)
                F_drift_y += du_dr * (dy / r)
                scalar_term += d2u_dr2 + (du_dr / r)
            end
        end
        E_kin += (F_drift_x^2 + F_drift_y^2 + scalar_term)
        F_drift_sq += F_drift_x^2 + F_drift_y^2
        Scalar_term_sum += scalar_term
    end
    return -0.5 * E_kin, 0.5*F_drift_sq,  -0.25 * Scalar_term_sum
end

"""
    local_energy(xcoord, ycoord, L, R_match, Constants) -> Tuple{Float64, Float64, Float64}

Calculates the local energy (kinetic + interaction) for a given configuration of particles. 

# Input:
- `xcoord::Vector{Float64}`: x-coordinates of all particles.
- `ycoord::Vector{Float64}`: y-coordinates of all particles.
- `L::Float64`: Box length (periodic boundary conditions assumed).
- `R_match::Float64`: Matching distance for the Jastrow factor.
- `Constants::Tuple{Float64, Float64, Float64}`: Tuple (C1, C2, C3) of Jastrow constants.

# Output:
- `Tuple{Float64, Float64, Float64, Float64}`: A tuple containing (E_total, E_kinetic, E_interaction).
"""
function local_energy(xcoord::Vector{Float64}, ycoord::Vector{Float64}, L::Float64, R_match::Float64, Constants::Tuple{Float64, Float64, Float64})::Tuple{Float64, Float64, Float64, Float64, Float64}
    E_kin, E_kin_drift, E_kin_laplacian = local_kinetic_energy(xcoord, ycoord, L, R_match, Constants)
    E_int = local_interaction_energy(xcoord, ycoord, L)
    return E_kin + E_int, # Total energy
           E_kin_drift + E_int, # Energy from drift estimator
           E_kin_laplacian + E_int, # Energy from laplacian estimator
           E_kin, # Kinetic energy only
           E_int # Interaction energy only
end

