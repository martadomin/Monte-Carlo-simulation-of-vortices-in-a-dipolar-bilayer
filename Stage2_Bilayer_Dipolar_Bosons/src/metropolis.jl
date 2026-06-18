include(normpath(joinpath(@__DIR__, "utils.jl")))
include(normpath(joinpath(@__DIR__, "jastrow.jl")))
include(normpath(joinpath(@__DIR__, "energy.jl")))

using Random, ProgressMeter, Plots

function move_one_part(x_A::Vector{Float64}, y_A::Vector{Float64},
                        x_B::Vector{Float64}, y_B::Vector{Float64},
                        delta::Float64, L::Float64)::Tuple{Int, Symbol, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}}

    x_A_new = copy(x_A)
    y_A_new = copy(y_A)
    x_B_new = copy(x_B)
    y_B_new = copy(y_B)

    if rand() < 0.5
        # Move a particle in layer A
        moved_id          = rand(1:length(x_A))
        x_A_new[moved_id] = wrap_position(x_A[moved_id] + rand() * 2 * delta - delta, L)
        y_A_new[moved_id] = wrap_position(y_A[moved_id] + rand() * 2 * delta - delta, L)
        return moved_id, :A, x_A_new, y_A_new, x_B_new, y_B_new
    else
        # Move a particle in layer B
        moved_id          = rand(1:length(x_B))
        x_B_new[moved_id] = wrap_position(x_B[moved_id] + rand() * 2 * delta - delta, L)
        y_B_new[moved_id] = wrap_position(y_B[moved_id] + rand() * 2 * delta - delta, L)
        return moved_id, :B, x_A_new, y_A_new, x_B_new, y_B_new
    end
end

function move_all_part(x_A, y_A, x_B, y_B, delta, L)
    x_A_new = wrap_position.(x_A .+ randn(length(x_A)) .* delta, L)
    y_A_new = wrap_position.(y_A .+ randn(length(y_A)) .* delta, L)
    x_B_new = wrap_position.(x_B .+ randn(length(x_B)) .* delta, L)
    y_B_new = wrap_position.(y_B .+ randn(length(y_B)) .* delta, L)
    return x_A_new, y_A_new, x_B_new, y_B_new
end

function compute_logΨ(x_A::Vector{Float64},
                      y_A::Vector{Float64},
                      x_B::Vector{Float64},
                      y_B::Vector{Float64},
                      R_match::Float64,
                      R0::Float64,
                      r_grid::Vector{Float64},
                      psi::Vector{Float64},
                      L::Float64,
                      Constants::Tuple{Float64, Float64, Float64})::Float64
    logΨ = 0.0
    N_half = length(x_A)

    # AA pairs
    @inbounds for i in 1:N_half
        @inbounds for j in (i + 1):N_half
            dx = get_periodic_difference(x_A[i], x_A[j], L)
            dy = get_periodic_difference(y_A[i], y_A[j], L)
            r = sqrt(dx^2 + dy^2)
            logΨ += u_AA(r, R_match, L, Constants)
        end
    end

    # BB pairs
    @inbounds for α in 1:N_half
        @inbounds for β in (α + 1):N_half
            dx = get_periodic_difference(x_B[α], x_B[β], L)
            dy = get_periodic_difference(y_B[α], y_B[β], L)
            r = sqrt(dx^2 + dy^2)
            logΨ += u_AA(r, R_match, L, Constants)
        end
    end

    # AB pairs
    @inbounds for i in 1:N_half
        @inbounds for α in 1:N_half
            dx = get_periodic_difference(x_A[i], x_B[α], L)
            dy = get_periodic_difference(y_A[i], y_B[α], L)
            r = sqrt(dx^2 + dy^2)
            logΨ += u_AB(r, R0, r_grid, psi)
        end
    end
    
    return logΨ
end

#Only valid when moving one particle for each step:
function compute_ΔlogΨ(x_A_old::Vector{Float64}, y_A_old::Vector{Float64},
                        x_B_old::Vector{Float64}, y_B_old::Vector{Float64},
                        x_A_new::Vector{Float64}, y_A_new::Vector{Float64},
                        x_B_new::Vector{Float64}, y_B_new::Vector{Float64},
                        id::Int, layer::Symbol,
                        R_match::Float64, L::Float64,
                        Constants::Tuple{Float64, Float64, Float64},
                        R0::Float64,
                        r_grid::Vector{Float64},
                        psi::Vector{Float64})::Float64

    ΔlogΨ  = 0.0
    N_half = length(x_A_old)

    if layer == :A
        @inbounds for j in 1:N_half

            # AA pairs: id vs other A particles
            if j != id
                r_old = sqrt(get_periodic_difference(x_A_old[id], x_A_old[j], L)^2 +
                             get_periodic_difference(y_A_old[id], y_A_old[j], L)^2)
                r_new = sqrt(get_periodic_difference(x_A_new[id], x_A_old[j], L)^2 +
                             get_periodic_difference(y_A_new[id], y_A_old[j], L)^2)
                ΔlogΨ += u_AA(r_new, R_match, L, Constants) -
                         u_AA(r_old, R_match, L, Constants)
            end

            # AB pairs: id (layer A) vs all B particles
            r_old = sqrt(get_periodic_difference(x_A_old[id], x_B_old[j], L)^2 +
                         get_periodic_difference(y_A_old[id], y_B_old[j], L)^2)
            r_new = sqrt(get_periodic_difference(x_A_new[id], x_B_old[j], L)^2 +
                         get_periodic_difference(y_A_new[id], y_B_old[j], L)^2)
            ΔlogΨ += u_AB(r_new, R0, r_grid, psi) -
                     u_AB(r_old, R0, r_grid, psi)
        end

    else  # layer == :B
        @inbounds for j in 1:N_half

            # BB pairs: id vs other B particles
            if j != id
                r_old = sqrt(get_periodic_difference(x_B_old[id], x_B_old[j], L)^2 +
                             get_periodic_difference(y_B_old[id], y_B_old[j], L)^2)
                r_new = sqrt(get_periodic_difference(x_B_new[id], x_B_old[j], L)^2 +
                             get_periodic_difference(y_B_new[id], y_B_old[j], L)^2)
                ΔlogΨ += u_AA(r_new, R_match, L, Constants) -
                         u_AA(r_old, R_match, L, Constants)
            end

            # AB pairs: id (layer B) vs all A particles
            r_old = sqrt(get_periodic_difference(x_B_old[id], x_A_old[j], L)^2 +
                         get_periodic_difference(y_B_old[id], y_A_old[j], L)^2)
            r_new = sqrt(get_periodic_difference(x_B_new[id], x_A_old[j], L)^2 +
                         get_periodic_difference(y_B_new[id], y_A_old[j], L)^2)
            ΔlogΨ += u_AB(r_new, R0, r_grid, psi) -
                     u_AB(r_old, R0, r_grid, psi)
        end
    end

    return ΔlogΨ
end

function tune_delta(
    x_A::Vector{Float64},
    y_A::Vector{Float64},
    x_B::Vector{Float64},
    y_B::Vector{Float64},
    L::Float64,
    R_match::Float64,
    Constants::Tuple{Float64, Float64, Float64},
    R0::Float64,
    r_grid::Vector{Float64},
    psi::Vector{Float64};
    target_ratio::Float64 = 0.5,
    num_tune_steps::Int = 10^4,
    block_size::Int = 100
)::Tuple{Float64, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}}

    delta = 0.1 * L 

    for _ in 1:num_tune_steps÷block_size
        accepted = 0
        for _ in 1:block_size
            moved_id, layer, x_A_new, y_A_new, x_B_new, y_B_new = move_one_part(x_A, y_A, x_B, y_B, delta, L)
            ΔlogΨ = compute_ΔlogΨ(x_A, y_A, x_B, y_B,x_A_new, y_A_new, x_B_new, y_B_new,moved_id, layer, R_match, L, Constants, R0, r_grid, psi)
            if log(rand()) < 2 * ΔlogΨ
                x_A = x_A_new
                y_A = y_A_new
                x_B = x_B_new
                y_B = y_B_new
                accepted += 1
            end
        end
        # adjust delta to push acceptance ratio towards target
        ratio = accepted / block_size
        delta *= (ratio + 1e-2) / target_ratio
        delta = min(delta, L/2)  
        delta = max(delta, 1e-6) 
    end

    return delta, x_A, y_A, x_B, y_B
end


"""
    metropolis(num_part, num_steps, delta, L, R_match, Constants, final_energy_plot, plot_every)

Runs a full Metropolis Monte Carlo simulation for the 2D dipolar system.

# Input:
- `num_part::Int`: Number of particles.
- `num_steps::Int`: Number of Metropolis steps.
- `delta::Float64`: Maximum displacement for particle moves.
- `L::Float64`: Length of the periodic simulation box.
- `R_match::Float64`: Matching radius used in the two-body Jastrow factor.
- `Constants::Tuple{Float64, Float64, Float64}`: Parameters used by the two-body correlation function.
- `final_energy_plot`: Controls whether the final energy plot is produced.
- `plot_every`: Interval used for plotting during the simulation.


# Notes
- Uses block averaging (`step_block`) for energy and structure factor sampling to reduce autocorrelation.
- Particle positions are stored and binned in [-L/2, L/2].
- Output: density and pair correlation histograms normalized as probability densities.
- The function displays a plot of the energy evolution over Monte Carlo steps.

# Usage
Call this function to simulate the equilibrium properties of the system, and to extract observables such as energy, density profiles, g2, and structure factor.

# Example
```julia
E_tot, E_sq, acceptance_ratio, E_kin, E_int, E_tot_drift = metropolis(30, 10^6, 0.1, sqrt(30), 0.5, (1.0, 1.0, 1.0); final_energy_plot=false, plot_every=100)
"""

function metropolis(
    num_part::Int,
    num_steps::Int,
    delta::Float64,
    L::Float64,
    h::Float64,
    R_match::Float64,
    R0::Float64,
    r_grid::Vector{Float64},
    psi::Vector{Float64},
    u_prime_grid::Vector{Float64},
    u_doubleprime_grid::Vector{Float64},
    Constants::Tuple{Float64, Float64, Float64};
    x_init::Union{Vector{Float64}, Nothing} = nothing,
    y_init::Union{Vector{Float64}, Nothing} = nothing,
    final_energy_plot::Bool = false,
    plot_every::Int = 100,
    progress::Bool = true,
    num_bins::Int = 100,
    move_all::Bool = false
    )::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Vector{Float64}, Vector{Float64}}

    acceptance_ratio = 0.0
    n_uncorr = 0
    step_block = 1
    E_tot = 0.0
    E_sq = 0.0
    E_kin = 0.0
    E_int = 0.0
    E_tot_drift = 0.0
    E_tot_laplacian = 0.0
    energies = Float64[]
    energies_drift = Float64[]
    energies_laplacian = Float64[]
    energy_plot = Float64[]
    step_plot = Int[]
    gr = zeros(Float64, num_bins)
    n_gr_samples = 0 
    E_local = NaN
    E_local_drift = NaN
    E_local_laplacian = NaN
    E_kinetic = NaN
    F_drift = NaN
    E_potential = NaN
    
    if x_init === nothing || y_init === nothing
        x_coord, y_coord = random_initial_config(num_part, L, "Uniform")
    else
        x_coord = copy(x_init)
        y_coord = copy(y_init)
    end

    x_A, y_A = x_coord[1:num_part÷2], y_coord[1:num_part÷2]
    x_B, y_B = x_coord[(num_part÷2 + 1):end], y_coord[(num_part÷2 + 1):end]

    _, _, E_local, E_local_drift, E_local_laplacian, E_kinetic, E_potential = energy_estimators(x_A, y_A, x_B, y_B, L, h, R_match, Constants, R0, r_grid, u_prime_grid, u_doubleprime_grid)
    @assert !isnan(E_local) "Initial configuration produced NaN energy. Re-initialize."
    
    if progress == true
        progress_bar = Progress(num_steps; desc="Running Metropolis $num_part...", showspeed=true)
    end

    logΨ_current = move_all ? compute_logΨ(x_A, y_A, x_B, y_B, R_match, R0, r_grid, psi, L, Constants) : 0.0

    for i in 1:num_steps
        if progress == true
            next!(progress_bar)
        end

        if move_all == true
            x_A_new, y_A_new, x_B_new, y_B_new = move_all_part(x_A, y_A, x_B, y_B, delta, L)
            logΨ_new = compute_logΨ(x_A_new, y_A_new, x_B_new, y_B_new, R_match, R0, r_grid, psi, L, Constants)
            ΔlogΨ = 2 * (logΨ_new - logΨ_current)
        else
            moved_id, layer, x_A_new, y_A_new, x_B_new, y_B_new = move_one_part(x_A, y_A, x_B, y_B, delta, L)
            ΔlogΨ = 2 * compute_ΔlogΨ(x_A, y_A, x_B, y_B, x_A_new, y_A_new, x_B_new, y_B_new, moved_id, layer, R_match, L, Constants, R0, r_grid, psi)
        end

        if log(rand()) < ΔlogΨ
            x_A = x_A_new
            y_A = y_A_new
            x_B = x_B_new
            y_B = y_B_new
            acceptance_ratio += 1

            if move_all
                logΨ_current = logΨ_new
            end

            _, _, E_local, E_local_drift, E_local_laplacian, E_kinetic, E_potential = energy_estimators(x_A, y_A, x_B, y_B, L, h, R_match, Constants, R0, r_grid, u_prime_grid, u_doubleprime_grid)
        end

        if i % step_block == 0

            if isnan(E_local)
                @warn "NaN detected at step $i: E_local = $E_local"
                continue
            end

            push!(energies, E_local)
            push!(energies_drift, E_local_drift)
            push!(energies_laplacian, E_local_laplacian)

            E_tot += E_local
            E_sq += E_local^2
            E_kin += E_kinetic
            E_int += E_potential
            E_tot_drift += E_local_drift
            E_tot_laplacian += E_local_laplacian

            n_uncorr += 1
            if i % plot_every == 0
                push!(step_plot, i)
                push!(energy_plot, E_local_drift / num_part)
            end
            x_coord = vcat(x_A, x_B)
            y_coord = vcat(y_A, y_B)
            n_gr_samples += 1
        end
    end

    x_coord = vcat(x_A, x_B)
    y_coord = vcat(y_A, y_B)

    if final_energy_plot
        p1 = plot(
        step_plot,
        energy_plot,
        xlabel="Step",
        ylabel="Local energy per particle",
        title="Energy evolution",
        legend=false,
        lw=2,
        )
        display(p1)
        println("\n>>> Energy Graph created. Press ENTER to continue...")
        readline()
    end

    println("Acceptance ratio: ", acceptance_ratio / num_steps)

    # r_vals, gr_normalized = normalize_gr(gr, num_part, L, n_gr_samples)

    return energies,
            energies_drift,
            energies_laplacian,
            E_tot / n_uncorr,
            E_sq / n_uncorr,
            E_tot_drift / n_uncorr,
            E_tot_laplacian / n_uncorr,
            E_kin / n_uncorr,
            E_int / n_uncorr,
            acceptance_ratio / num_steps,
            x_coord,
            y_coord
end