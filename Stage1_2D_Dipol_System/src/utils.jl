using Random

"""
    get_periodic_difference(x1::Float64, x2::Float64, L::Float64) -> Float64

Computes the minimum-image (periodic) difference between two points in a 1D periodic box.

# Input:
- `x1::Float64`: Position of the first point.
- `x2::Float64`: Position of the second point.
- `L::Float64`: Length of the periodic box.

# Output:
- `Float64`: The difference `(x1 - x2)`, mapped to the interval [0, L].

# Notes
Useful for applying periodic boundary conditions and minimum-image convention in simulations.
"""
function get_periodic_difference(x1::Float64, x2::Float64, L::Float64)::Float64
    diff = x1 - x2
    # Shift to [0, L), then to [-L/2, L/2]
    return diff - L * round(diff / L)
end

"""
    random_initial_config(num_part::Int, L::Float64) -> Matrix{Float64}

Generates a random initial configuration of `num_part` particles uniformly distributed in a 1D periodic box of length `L`.

# Input:
- `num_part::Int`: Number of particles.
- `L::Float64`: Length of the simulation box.

# Output:
- `Matrix{Float64}`: Positions of all particles, each in the interval [-L/2, L/2].

# Notes
Particle positions are initialized randomly and independently with uniform probability over the full simulation box.
"""
function random_initial_config(num_part::Int, L::Float64, distribution::AbstractString)::Tuple{Vector{Float64}, Vector{Float64}}
    if distribution == "Uniform"
        positions = L .* rand(2, num_part) .- L/2  # Random initial configuration in the range [-L/2, L/2]
    else
        error("Unsupported distribution: $distribution")
    end

    return positions[1, :], positions[2, :]  # Return x and y coordinates as separate vectors

end
