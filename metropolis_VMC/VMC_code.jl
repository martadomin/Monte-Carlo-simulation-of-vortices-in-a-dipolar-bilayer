using Random, LinearAlgebra, Dierckx, StatsBase, Plots, Base.Threads, LaTeXStrings, ProgressMeter, Roots, Bessels


# ------------------------------------------
# Utility Functions
# ------------------------------------------
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
    if r < 1e-10
        return 0.0  # or some safe fallback
    end

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
    if r < 1e-10
        return 0.0
    end
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
    if r < 1e-10
        return 0.0
    end

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

# ========== Energy Calculations ==========

"""
    interaction_energy(x_coord, num_part, V0, k_lat) -> Float64

Calculates the total interaction energy for a set of particles, each pair interacting via a cavity-mediated cosine potential.

# Input:
- `positions::Matrix{Float64}`: Positions of all particles.
- `L::Float64`: Length of the periodic box.

# Output:
- `Float64`: The total potential energy, sum over all unique pairs.

# Notes
The interaction is given by (|r_i - r_j|^(-3)) for each pair (i < j).
"""
function local_interaction_energy(x_coord::Vector{Float64}, y_coord::Vector{Float64}, L::Float64)::Float64
    E_int = 0.0
    num_part = length(x_coord)

    @inbounds for i in 1:num_part
        @inbounds for j in (i + 1):num_part
            dx = get_periodic_difference(x_coord[i], x_coord[j], L)
            dy = get_periodic_difference(y_coord[i], y_coord[j], L)
            E_int += (dx^2 + dy^2)^(-3/2)
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
                                Constants::Tuple{Float64, Float64, Float64})::Float64
    
    num_part = length(xcoord)
    E_kin = 0.0

    @inbounds for k in 1:num_part
        F_drift_x = 0.0
        F_drift_y = 0.0
        scalar_term = 0.0
        for i in 1:num_part
            if i != k
                dx = get_periodic_difference(xcoord[i], xcoord[k], L)
                dy = get_periodic_difference(ycoord[i], ycoord[k], L)
                r = sqrt(dx^2 + dy^2)
                du_dr = u2_first_derivative(r, R_match, L, Constants)
                d2u_dr2 = u2_second_derivative(r, R_match, L, Constants)
                F_drift_x += du_dr * (dx/r)
                F_drift_y += du_dr * (dy/r)
                scalar_term += d2u_dr2 + (du_dr/r)
            end
        end
        E_kin += (F_drift_x^2 + F_drift_y^2 + scalar_term)
    end

    return -0.5 * E_kin
end

function local_energy(xcoord::Vector{Float64}, ycoord::Vector{Float64}, L::Float64, R_match::Float64, Constants::Tuple{Float64, Float64, Float64})::Tuple{Float64, Float64, Float64}
    E_kin = local_kinetic_energy(xcoord, ycoord, L, R_match, Constants)
    E_int = local_interaction_energy(xcoord, ycoord, L)
    return E_kin + E_int, E_kin, E_int
end

# ------------------------------------------
# Creating the Wavefunction
# ------------------------------------------

function compute_logΨ(xcoord::Vector{Float64}, ycoord::Vector{Float64}, R_match::Float64, L::Float64, Constants::Tuple{Float64, Float64, Float64})::Float64
    logΨ = 0.0
    num_part = length(xcoord)

    @inbounds for i in 1:num_part
        @inbounds for j in (i + 1):num_part
            dx = get_periodic_difference(xcoord[i], xcoord[j], L)
            dy = get_periodic_difference(ycoord[i], ycoord[j], L)
            r = sqrt(dx^2 + dy^2)
            logΨ += u2(r, R_match, L, Constants)
        end
    end
    return logΨ
end

#Only valid when moving one particle for each step:
function compute_ΔlogΨ(x_old::Vector{Float64}, y_old::Vector{Float64}, x_new::Vector{Float64}, y_new::Vector{Float64}, R_match::Float64, L::Float64, Constants::Tuple{Float64, Float64, Float64}, id::Int)::Float64
    ΔlogΨ = 0.0
    num_part = length(x_old)

    @inbounds for i in 1:num_part
        if i != id
            #old distance
            dx_old = get_periodic_difference(x_old[i], x_old[id], L)
            dy_old = get_periodic_difference(y_old[i], y_old[id], L)
            r_old = sqrt(dx_old^2 + dy_old^2)
            #new distance
            dx_new = get_periodic_difference(x_new[i], x_new[id], L)
            dy_new = get_periodic_difference(y_new[i], y_new[id], L)
            r_new = sqrt(dx_new^2 + dy_new^2)
            
            ΔlogΨ += u2(r_new, R_match, L, Constants) - 
                        u2(r_old, R_match, L, Constants)
        end
    end
    return ΔlogΨ
end

# ------------------------------------------
# System Initialization and Particle Moves
# ------------------------------------------

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
function random_initial_config(num_part::Int, L::Float64)::Tuple{Vector{Float64}, Vector{Float64}}
    positions = L .* rand(2, num_part) .- L/2  # Random initial configuration in the range [-L/2, L/2]
    return positions[1, :], positions[2, :]  # Return x and y coordinates as separate vectors
end

"""
    move_one_part(x_coord::Vector{Float64}, y_coord::Vector{Float64}, num_part::Int, delta::Float64, L::Float64) -> Tuple{Int, Vector{Float64}, Vector{Float64}}

Proposes a random move for a single particle by displacing it within [-delta, delta] and applying periodic boundary conditions.

# Input:
- `x_coord::Vector{Float64}`: Current positions of all particles in the x direction.
- `y_coord::Vector{Float64}`: Current positions of all particles in the y direction.
- `num_part::Int`: Number of particles.
- `delta::Float64`: Maximum displacement (half-width of the move interval).
- `L::Float64`: Length of the simulation box.

# Output:
- `Int`: Index of the particle that was moved.
- `Vector{Float64}`: New positions in the x direction after the proposed move (with periodicity).
- `Vector{Float64}`: New positions in the y direction after the proposed move (with periodicity).

# Notes
A particle is selected at random and displaced by a random amount in [-delta, delta]. The new position is wrapped to the periodic box using the minimum image convention.
"""
function move_one_part(x_coord::Vector{Float64}, y_coord::Vector{Float64}, delta::Float64, L::Float64)::Tuple{Int, Vector{Float64}, Vector{Float64}}
    x_coord_new = copy(x_coord)
    y_coord_new = copy(y_coord)
    id = rand(1:length(x_coord))
    x_coord_new[id] += rand() * (2 * delta) - delta  # Displacement in [-delta, delta]
    x_coord_new[id] = get_periodic_difference(x_coord_new[id], 0.0, L)  # Apply periodic boundary conditions
    y_coord_new[id] += rand() * (2 * delta) - delta  # Displacement in [-delta, delta]
    y_coord_new[id] = get_periodic_difference(y_coord_new[id], 0.0, L)  # Apply periodic boundary conditions
    return id, x_coord_new, y_coord_new
end

# ------------------------------------------
# Final Metropolis Implementation
# ------------------------------------------

"""
    metropolis(num_part, num_steps, num_bins, delta, L, V0, k_lat, psi_interp, k_L, k_contact, α, long_range, fermi_stats, reatto_chester, contact)
        -> Tuple{Float64, Float64, Vector{ComplexF64}, Vector{Float64}, Matrix{Float64}, Float64, Vector{Float64}}

Runs a full Metropolis Monte Carlo simulation for a 1D quantum system with customizable two-body wavefunction structure and interactions.

# Input:
- `num_part::Int`: Number of particles.
- `num_steps::Int`: Number of Metropolis steps.
- `num_bins::Int`: Number of bins for histograms (density, pair).
- `delta::Float64`: Maximum displacement for particle moves.
- `L::Float64`: Length of the periodic simulation box.

# Output:
- `Float64`: Mean total energy per configuration (E_tot / n_uncorr).
- `Float64`: Mean squared total energy per configuration.
- `Vector{ComplexF64}`: Static structure factor S(k) (for a range of k values).
- `Vector{Float64}`: Normalized 1D density histogram (n(x)).
- `Matrix{Float64}`: Normalized 2D pair density histogram (g2(x, x')).
- `Float64`: Acceptance ratio of proposed moves.
- `Vector{Float64}`: Block-averaged local energy per particle, for convergence diagnostics.

# Notes
- Uses block averaging (`step_block`) for energy and structure factor sampling to reduce autocorrelation.
- Particle positions are stored and binned in [-L/2, L/2].
- Output: density and pair correlation histograms normalized as probability densities.
- The function displays a plot of the energy evolution over Monte Carlo steps.

# Usage
Call this function to simulate the equilibrium properties of the system, and to extract observables such as energy, density profiles, g2, and structure factor.

# Example
```julia
E, E2, SSF, n_x, g2_xx, acc_ratio, E_trace = metropolis(8, 10^6, 100, 0.05, 1.0, 1.0, 2π, psi_interp, 2.0, 1.0, 1.0, true, true, false, false)
"""
function metropolis(num_part::Int, num_steps::Int, num_bins::Int, delta::Float64, L::Float64, R_match::Float64, Constants::Tuple{Float64, Float64, Float64}; live_plot::Bool=false, plot_every::Int=100)
    acceptance_ratio = 0.0
    n_uncorr = 0
    # bins = range(-L/2, stop=L/2, length=num_bins+1)
    # hist_1d = zeros(Float64, num_bins)
    # hist_2d = zeros(Float64, num_bins, num_bins)
    # dx = L / num_bins
    step_block = 1
    E_tot = 0.0
    E_sq = 0.0
    energy_trace = Float64[]
    step_trace = Int[]
    # E_local_values = Vector{Float64}(undef, (num_steps÷step_block))
    # iter_val = Vector{Float64}(undef, (num_steps÷step_block))
    # idx_plot = 1
    # configurations = Vector{Vector{Float64}}(undef, num_steps ÷ step_block)

    # final_point = 5*L
    # k = (2*π/L) * collect(1:1:final_point)
    # SSF = zeros(ComplexF64, length(k))

    x_coord, y_coord = random_initial_config(num_part, L)
    
    progress = Progress(num_steps; desc="Running Metropolis $num_part...", showspeed=true)
    for i in 1:num_steps
        next!(progress)  # Update progress bar

        moved_id, x_new, y_new = move_one_part(x_coord, y_coord, delta, L)

        ΔlogΨ = compute_ΔlogΨ(x_coord, y_coord, x_new, y_new, R_match, L, Constants, moved_id)

        if log(rand()) < 2 * ΔlogΨ
            x_coord = x_new
            y_coord = y_new
            acceptance_ratio += 1
        end
    
        # for x in x_coord
        #     bin_idx = min(num_bins, max(1, Int(floor((x + L/2) / L * num_bins)) + 1))
        #     hist_1d[bin_idx] += 1
        # end
    
        # for i in 1:num_part
        #     for j in (i + 1):num_part
        #         bin_x = min(num_bins, max(1, Int(floor((x_coord[i] + L/2) / L * num_bins)) + 1))
        #         bin_y = min(num_bins, max(1, Int(floor((x_coord[j] + L/2) / L * num_bins)) + 1))

        #         hist_2d[bin_x, bin_y] += 1
        #     end
        # end
    
        if i % step_block == 0

            E_local, E_kinetic, E_potential = local_energy(x_coord, y_coord, L, R_match, Constants)
            
            if isnan(E_local)
                @warn "NaN detected at step $i: E_local = $E_local"
                continue  # Skip this iteration to avoid polluting data
            end

            E_tot += E_local
            E_sq += E_local^2
            n_uncorr += 1
            if i % plot_every == 0
                push!(step_trace, i)
                push!(energy_trace, E_local / num_part)

            end
            
            # for a in 1:num_part, b in 1:num_part
            #     SSF .+= exp.(im * (x_coord[a] - x_coord[b]) .* k)
            # end
    
            # iter_val[idx_plot] = i
            # E_local_values[idx_plot] = E_local / num_part
            # configurations[idx_plot] = copy(x_coord)
            # idx_plot += 1
        end
    end

    if live_plot==true
        energy_plot = plot(
        step_trace,
        energy_trace,
        xlabel="Step",
        ylabel="Local energy per particle",
        title="Energy evolution",
        legend=false,
        lw=2,
        )
        display(energy_plot)
        readline()
    end

    println("Acceptance ratio: ", acceptance_ratio / num_steps)
    # hist_1d ./= (sum(hist_1d) * dx)
    # hist_2d ./= (sum(hist_2d) * dx^2)
    # SSF ./= n_uncorr

    return E_tot / n_uncorr, E_sq / n_uncorr, acceptance_ratio / num_steps
end


#System parameters
n = 32
num_part = 5
L = sqrt(num_part/n)

println("Density n = ", n)
println("Number of particles: ", num_part)
println("L = ", L)

#Simulation parameters
R_match_vals = LinRange(0.01, L, 100) # Matching radius to validate continuity of f2

results_path = joinpath(@__DIR__, "R_match_energy_results.txt")
R_match_vals_vec = collect(R_match_vals)
energy_means = Vector{Float64}(undef, length(R_match_vals_vec))
energy_sq_means = Vector{Float64}(undef, length(R_match_vals_vec))
acceptance_ratio_vec = Vector{Float64}(undef, length(R_match_vals_vec))

@threads for idx in eachindex(R_match_vals_vec)
    R_match = R_match_vals_vec[idx]
    num_steps = 10^6
    num_bins = 100
    delta = 0.1
    Constants = calculate_constants(L, R_match)
    E_mean, E_sq_mean, acceptance_ratio = metropolis(num_part, num_steps, num_bins, delta, L, R_match, Constants; live_plot=false, plot_every=10^3)
    energy_means[idx] = E_mean
    energy_sq_means[idx] = E_sq_mean
    acceptance_ratio_vec[idx] = acceptance_ratio
end

open(results_path, "w") do io
    println(io, "R_match\tMeanEnergyPerParticle\tMeanSquaredEnergy\tAcceptanceRatio")
    for idx in eachindex(R_match_vals_vec)
        println(io, "$(R_match_vals_vec[idx])\t$(energy_means[idx])\t$(energy_sq_means[idx])\t$(acceptance_ratio_vec[idx])")
    end
end

println("Saved sweep results to: ", results_path)
