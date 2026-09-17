using Random

"""
    get_periodic_difference(x1::Float64, x2::Float64, L::Float64) -> Float64
Computes the minimum-image (periodic) difference between two points in a 1D periodic box.
# Input:
- `x1::Float64`: Position of the first point.
- `x2::Float64`: Position of the second point.
- `L::Float64`: Length of the periodic box.
# Output:
- `Float64`: The difference `(x1 - x2)`, mapped to the interval [-L/2, L/2).
"""
function get_periodic_difference(x1::Float64, x2::Float64, L::Float64)::Float64
    return mod(x1 - x2 + L/2, L) - L/2
end

"""
    wrap_position(x::Float64, L::Float64) -> Float64
Wraps a position back into the simulation box [0, L).
# Input:
- `x::Float64`: Position to wrap.
- `L::Float64`: Length of the periodic box.
# Output:
- `Float64`: Position mapped to [0, L).
"""
function wrap_position(x::Float64, L::Float64)::Float64
    return mod(x, L)
end

"""
    random_initial_config(num_part::Int, L::Float64, distribution::AbstractString) -> Tuple
Generates a random initial configuration of `num_part` particles in a 2D periodic box.
# Input:
- `num_part::Int`: Number of particles.
- `L::Float64`: Length of the simulation box.
- `distribution::AbstractString`: Type of distribution ("Uniform").
# Output:
- `Tuple{Vector{Float64}, Vector{Float64}}`: x and y coordinates in [0, L).
"""
function random_initial_config(num_part::Int, L::Float64, distribution::AbstractString)::Tuple{Vector{Float64}, Vector{Float64}}
    if distribution == "Uniform"
        positions = L .* rand(2, num_part)  # [0, L)
    elseif distribution == "Normal"
        positions = L/2 .+ (L/4) .* randn(2, num_part)  # Centered at L/2 with std dev L/4
        positions = wrap_position.(positions, L)  # Wrap into [0, L)    
    else
        error("Unsupported distribution: $distribution")
    end
    return positions[1, :], positions[2, :]
end