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
                if r_ij > 1e-10
                    du_dr = u2_first_derivative(r_ij, R_match, L, Constants)
                    drift_x +=  du_dr * (dx / r_ij)
                    drift_y +=  du_dr * (dy / r_ij)
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

function weight_update(E_loc_old::Float64, E_loc_new::Float64, E_ref::Float64, Δτ::Float64)::Float64
    return exp(- Δτ * (0.5*(E_loc_old + E_loc_new) - E_ref))
end

function branching_step(weights::Vector{Float64})::Vector{Int}
    return [floor(Int, weight + rand()) for weight in weights]
end

function population_control(num_walkers::Int, num_target::Int, average_gs_E::Float64, Δτ::Float64)::Float64
    return average_gs_E - (1 / Δτ) * log(num_walkers / num_target)
end

function dmc(x_init::Vector{Float64}, y_init::Vector{Float64},
             num_walkers::Int, num_part::Int, num_steps::Int,
             Δτ::Float64, L::Float64, R_match::Float64,
             Constants::Tuple{Float64, Float64, Float64},
             E_ref_initial::Float64, num_target::Int;
             num_equil::Int = num_steps ÷ 5,
             plot_energy::Bool = false)

    D = 0.5  # dimensionless units

    # Initialize walkers from VMC config
    x_walkers = [copy(x_init) for _ in 1:num_walkers]
    y_walkers = [copy(y_init) for _ in 1:num_walkers]
    weights   = ones(Float64, num_walkers)
    E_ref     = E_ref_initial
    E_history = Float64[]
    E_plot = Float64[]

    E_loc_old = Vector{Float64}(undef, num_walkers)
    @threads for i in 1:length(x_walkers)
        E_loc_old[i], _, _, _, _ = local_energy(x_walkers[i], y_walkers[i], L, R_match, Constants)
    end

    prog = Progress(num_steps; desc="Running DMC...", showspeed=true)
    for step in 1:num_steps

        new_x = Vector{Vector{Float64}}(undef, length(x_walkers))
        new_y = Vector{Vector{Float64}}(undef, length(y_walkers))
        next!(prog)
        # 1. Drift Δt/2 at old positions
        @threads for i in 1:length(x_walkers)
            new_x[i], new_y[i] = drift_step(x_walkers[i], y_walkers[i], Δτ/2, L, R_match, Constants)
        end

        # 2. Diffusion Δt at drifted positions
        @threads for i in 1:length(new_x)
            new_x[i], new_y[i] = diffusion_step(new_x[i], new_y[i], L, D, Δτ)
        end

        # 3. Drift Δt/2 at new positions
        @threads for i in 1:length(x_walkers)
            new_x[i], new_y[i] = drift_step(new_x[i], new_y[i], Δτ/2, L, R_match, Constants)
        end

        x_walkers = new_x
        y_walkers = new_y

        # 4. Local energies and weights
        E_loc_new = Vector{Float64}(undef, length(x_walkers))
        @threads for i in 1:length(x_walkers)
            E_loc_new[i], _, _, _, _ = local_energy(x_walkers[i], y_walkers[i], L, R_match, Constants)
        end

        # 5. Branching weights update
        for i in 1:length(x_walkers)
            weights[i] *= weight_update(E_loc_old[i], E_loc_new[i], E_ref, Δτ)
        end

        avg_E = mean(E_loc_new)

        # 6. Branching
        num_copies = branching_step(weights)

        num_copies = branching_step(weights)
        new_x, new_y = Vector{Vector{Float64}}(), Vector{Vector{Float64}}()
        E_loc_old = Float64[]
        for i in 1:length(x_walkers)
            for _ in 1:num_copies[i]
                push!(new_x, copy(x_walkers[i]))
                push!(new_y, copy(y_walkers[i]))
                push!(E_loc_old, E_loc_new[i])  # carry forward with walker
            end
        end

        x_walkers = new_x
        y_walkers = new_y
        weights   = ones(Float64, length(x_walkers))

        if isempty(x_walkers)
            @warn "All walkers died at step $step"
            break
        end

        # 7. Population control
        E_ref = population_control(length(x_walkers), num_target, avg_E, Δτ)

        # 8. Accumulate after equilibration
        if step > num_equil
            push!(E_history, avg_E)
        end

        E_plot = push!(E_plot, avg_E)

        if step % 10000 == 0
            println("\nCheckpoint at step $step: Walkers = $(length(x_walkers)), Avg E_loc = $avg_E, E_ref = $E_ref")
            println("num_copies stats: min=$(minimum(num_copies)), max=$(maximum(num_copies)), mean=$(mean(num_copies))")
        end

    end

    if plot_energy
        p_energy = plot(1:length(E_plot), E_plot,
             xlabel="DMC step", ylabel="E", title="DMC Energy History", linewidth=2)
        display(p_energy)
        println("\n>>> DMC Energy History plot created. Press ENTER to continue...")
        readline()
    end

    E_dmc     = mean(E_history)
    E_dmc_err = std(E_history) / sqrt(length(E_history))
    
    return E_dmc, E_dmc_err, E_history
end
