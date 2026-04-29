function diffusion_step(x_coord::Vector{Float64}, y_coord::Vector{Float64}, L::Float64, D::Float64, Δτ::Float64)::Tuple{Vector{Float64}, Vector{Float64}}
    # Propose new positions with diffusion
    x_proposed = x_coord .+ sqrt(2 * D * Δτ) .* randn(length(x_coord))
    y_proposed = y_coord .+ sqrt(2 * D * Δτ) .* randn(length(y_coord))

    # Wrap proposed positions back into the box
    x_proposed = wrap_position.(x_proposed, L)
    y_proposed = wrap_position.(y_proposed, L)
    
    return x_proposed, y_proposed
end

function drift_step(x_coord::Vector{Float64}, y_coord::Vector{Float64}, Δτ::Float64, L::Float64, R_match::Float64, Constants::Tuple{Float64, Float64, Float64})::Tuple{Vector{Float64}, Vector{Float64}}
    num_part = length(x_coord)
    x_drifted = copy(x_coord)
    y_drifted = copy(y_coord)
    
    for i in 1:num_part
        drift_x = 0.0
        drift_y = 0.0
        for j in 1:num_part
            if i != j
                dx = get_periodic_difference(x_coord[i], x_coord[j], L)
                dy = get_periodic_difference(y_coord[i], y_coord[j], L)
                r_ij = sqrt(dx^2 + dy^2)
                if r_ij < L/2
                    du_dr = u2_first_derivative(r_ij, R_match, L, Constants)
                    if r_ij > 1e-10
                        drift_x += du_dr * (dx / r_ij)
                        drift_y += du_dr * (dy / r_ij)
                    end
                end
            end
        end
        x_drifted[i] += drift_x * Δτ
        y_drifted[i] += drift_y * Δτ
        
        # Wrap drifted positions back into the box
        x_drifted[i] = wrap_position(x_drifted[i], L)
        y_drifted[i] = wrap_position(y_drifted[i], L)
    end
    
    return x_drifted, y_drifted
end

function weight_update(E_loc::Float64, E_ref::Float64, Δτ::Float64)::Float64
    return exp(- Δτ * (E_loc - E_ref))
end

function branching_step(weights::Vector{Float64})::Vector{Int}
    return [floor(Int, weight + rand()) for weight in weights]
end

function population_control(num_walkers::Int, num_target::Int, average_gs_E::Float64, Δτ::Float64)::Float64
    return average_gs_E - (1 / Δτ) * log(num_walkers / num_target)
end

function dmc(num_walkers::Int, num_steps::Int, Δτ::Float64, L::Float64, R_match::Float64, Constants::Tuple{Float64, Float64, Float64}, E_ref_initial::Float64)
    # Initialize walkers
    x_walkers = [L .* rand(length(Constants)) for _ in 1:num_walkers]
    y_walkers = [L .* rand(length(Constants)) for _ in 1:num_walkers]
    weights = ones(num_walkers)
    E_ref = E_ref_initial
    E = 0.0

    for step in 1:num_steps
        # Propose new positions with diffusion and drift
        for i in 1:num_walkers
            x_diffused, y_diffused = diffusion_step(x_walkers[i], y_walkers[i], L, D, Δτ)
            x_drifted, y_drifted = drift_step(x_diffused, y_diffused, Δτ, L, R_match, Constants)
            x_walkers[i] = x_drifted
            y_walkers[i] = y_drifted
        end
        
        # Calculate local energies and update weights
        for i in 1:num_walkers
            E_loc = local_energy(x_walkers[i], y_walkers[i], L, R_match, Constants)
            E += E_loc
            weights[i] *= weight_update(E_loc, E_ref, Δτ)
        end
        
        # Branching step
        num_copies = branching_step(weights)
        
        # Create new walker population based on branching
        new_x_walkers = Vector{Vector{Float64}}()
        new_y_walkers = Vector{Vector{Float64}}()
        for i in 1:num_walkers
            for _ in 1:num_copies[i]
                push!(new_x_walkers, copy(x_walkers[i]))
                push!(new_y_walkers, copy(y_walkers[i]))
            end
        end
        
        x_walkers = new_x_walkers
        y_walkers = new_y_walkers
        
        # Population control
        if length(x_walkers) > 0
            E_ref = population_control(num_walkers, num_target=1000, average_gs_E=E/step, Δτ)
        else
            println("All walkers died at step $step")
            break
        end