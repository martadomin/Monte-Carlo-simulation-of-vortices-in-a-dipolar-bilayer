using Random

"""
    get_periodic_difference(x1::Float64, x2::Float64, L::Float64) -> Float64

Computes the minimum-image (periodic) difference between two points in a 1D periodic box.

# Input:
- `x1::Float64`: Position of the first point.
- `x2::Float64`: Position of the second point.
- `L::Float64`: Length of the periodic box.

# Output:
- `Float64`: The difference `(x1 - x2)`, mapped to the interval [-L/2, L/2].

# Notes
Useful for applying periodic boundary conditions and minimum-image convention in simulations.
"""
function get_periodic_difference(x1::Float64, x2::Float64, L::Float64)::Float64
    diff = x1 - x2
    # Shift to [0, L), then to [-L/2, L/2]
    return mod(diff + L/2, L) - L/2
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

function tune_delta(
    x_coord::Vector{Float64},
    y_coord::Vector{Float64},
    L::Float64,
    R_match::Float64,
    Constants::Tuple{Float64, Float64, Float64};
    target_ratio::Float64 = 0.5,
    num_tune_steps::Int = 10^4,
    block_size::Int = 100
)::Tuple{Float64, Vector{Float64}, Vector{Float64}}

    delta = 0.1 * L  # initial guess

    for _ in 1:num_tune_steps÷block_size
        accepted = 0
        for _ in 1:block_size
            moved_id, x_new, y_new = move_one_part(x_coord, y_coord, delta, L)
            ΔlogΨ = compute_ΔlogΨ(x_coord, y_coord, x_new, y_new,
                                   R_match, L, Constants, moved_id)
            if log(rand()) < 2 * ΔlogΨ
                x_coord = x_new
                y_coord = y_new
                accepted += 1
            end
        end
        # adjust delta to push acceptance ratio towards target
        ratio = accepted / block_size
        delta *= ratio / target_ratio
        delta = min(delta, L/2)  
        delta = max(delta, 1e-6) 
    end

    println("Tuned delta = $delta")
    return delta, x_coord, y_coord
end