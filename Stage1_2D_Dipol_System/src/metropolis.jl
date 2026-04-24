include(normpath(joinpath(@__DIR__, "utils.jl")))
include(normpath(joinpath(@__DIR__, "jastrow.jl")))
include(normpath(joinpath(@__DIR__, "energy.jl")))

using Random, ProgressMeter, Plots

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

    return delta, x_coord, y_coord
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
    R_match::Float64,
    Constants::Tuple{Float64, Float64, Float64};
    x_init::Union{Vector{Float64}, Nothing} = nothing,
    y_init::Union{Vector{Float64}, Nothing} = nothing,
    final_energy_plot::Bool = false,
    plot_every::Int = 100,
    progress::Bool = true
    )::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Float64, Float64, Float64, Float64, Float64, Float64, Float64}

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
    energy_plot_1 = Float64[]
    energy_plot_2 = Float64[]
    step_plot = Int[]
    
    if x_init === nothing || y_init === nothing
        x_coord, y_coord = random_initial_config(num_part, L, "Uniform")
    else
        x_coord = copy(x_init)
        y_coord = copy(y_init)
    end
    
    if progress == true
        progress_bar = Progress(num_steps; desc="Running Metropolis $num_part...", showspeed=true)
    end

    for i in 1:num_steps
        if progress == true
            next!(progress_bar)
        end

        moved_id, x_new, y_new = move_one_part(x_coord, y_coord, delta, L)

        ΔlogΨ = compute_ΔlogΨ(x_coord, y_coord, x_new, y_new, R_match, L, Constants, moved_id)

        if log(rand()) < 2 * ΔlogΨ
            x_coord = x_new
            y_coord = y_new
            acceptance_ratio += 1
        end
        
        if i % step_block == 0
            E_local, E_local_drift, E_local_laplacian, E_kinetic, E_potential = local_energy(x_coord, y_coord, L, R_match, Constants)
            push!(energies, E_local)
            push!(energies_drift, E_local_drift)
            push!(energies_laplacian, E_local_laplacian)

            if isnan(E_local)
                @warn "NaN detected at step $i: E_local = $E_local"
                continue  # Skip this iteration to avoid polluting data
            end

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
                push!(energy_plot_1, E_local_laplacian / num_part)
                push!(energy_plot_2, E_local / num_part)

            end
        end
    end

    if final_energy_plot
        energy_plot = plot(
        step_plot,
        energy_plot,
        xlabel="Step",
        ylabel="Local energy per particle",
        title="Energy evolution",
        legend=false,
        lw=2,
        )
        display(energy_plot)
        # Si ejecutas desde terminal, esto evita que el script termine y cierre el plot
        println("\n>>> Gráfico de energía generado. Presiona ENTER para continuar...")
        if !isinteractive() # Solo bloquea si no estás en un REPL interactivo
            readline()
        end
    end

    println("Acceptance ratio: ", acceptance_ratio / num_steps)

    return energies,
            energies_drift,
            energies_laplacian,
            E_tot / n_uncorr,
            E_sq / n_uncorr,
            E_tot_drift / n_uncorr,
            E_tot_laplacian / n_uncorr,
            E_kin / n_uncorr,
            E_int / n_uncorr,
            acceptance_ratio / num_steps
end