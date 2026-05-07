function diffusion_step(x_coord::Vector{Float64}, y_coord::Vector{Float64}, L::Float64, D::Float64, Δτ::Float64)::Tuple{Vector{Float64}, Vector{Float64}}
    # Propose new positions with diffusion
    x_proposed = x_coord .+ sqrt(2 * D * Δτ) .* randn(length(x_coord))
    y_proposed = y_coord .+ sqrt(2 * D * Δτ) .* randn(length(y_coord))

    # Wrap proposed positions back into the box
    x_proposed = wrap_position.(x_proposed, L)
    y_proposed = wrap_position.(y_proposed, L)
    
    return x_proposed, y_proposed
end

# function drift_step(x_coord::Vector{Float64}, y_coord::Vector{Float64}, Δτ::Float64, L::Float64, R_match::Float64, Constants::Tuple{Float64, Float64, Float64})::Tuple{Vector{Float64}, Vector{Float64}}
#     num_part = length(x_coord)
#     x_drifted = copy(x_coord)
#     y_drifted = copy(y_coord)
    
#     for i in 1:num_part
#         drift_x = 0.0
#         drift_y = 0.0
#         for j in 1:num_part
#             if i != j
#                 dx = get_periodic_difference(x_coord[i], x_coord[j], L)
#                 dy = get_periodic_difference(y_coord[i], y_coord[j], L)
#                 r_ij = sqrt(dx^2 + dy^2)
#                 if r_ij > 1e-10
#                     du_dr = u2_first_derivative(r_ij, R_match, L, Constants)
#                     drift_x +=  du_dr * (dx / r_ij)
#                     drift_y +=  du_dr * (dy / r_ij)
#                 end
#             end
#         end

#         x_drifted[i] += drift_x * Δτ
#         y_drifted[i] += drift_y * Δτ
        
#         # Wrap drifted positions back into the box
#         x_drifted[i] = wrap_position(x_drifted[i], L)
#         y_drifted[i] = wrap_position(y_drifted[i], L)
#     end
    
#     return x_drifted, y_drifted
# end

function weight_update(E_loc_old::Float64, E_loc_new::Float64, E_ref::Float64, Δτ::Float64)::Float64

    if isnan(E_loc_old) || isnan(E_loc_new)
        return 0.0 
    end

    exponent = - Δτ * (0.5*(E_loc_old + E_loc_new) - E_ref)
    exponent_capped = clamp(exponent, -50.0, log(10.0))
    
    return exp(exponent_capped)
end

function branching_step(weights::Vector{Float64})::Vector{Int}
    # Catch NaNs, ensure non-negative, and strictly limit max copies to 3
    safe_weights = [isnan(w) ? 0.0 : clamp(w, 0.0, 3.0) for w in weights]
    
    return [floor(Int, weight + rand()) for weight in safe_weights]
end

function population_control(num_walkers::Int, num_target::Int, avg_E_loc::Float64, Δτ::Float64)::Float64
    return avg_E_loc - (1 / Δτ) * log(num_walkers / num_target)
end

function dmc(x_init::Vector{Float64}, y_init::Vector{Float64},
             num_walkers::Int, num_part::Int, num_steps::Int,
             Δτ::Float64, L::Float64, R_match::Float64,
             Constants::Tuple{Float64, Float64, Float64},
             E_ref_initial::Float64, num_target::Int;
             num_equil::Int = num_steps ÷ 5,
             quadratic::Bool = false,
             plot_energy::Bool = false)

    D = 0.5

    # Initialize walkers
    x_walkers = [copy(x_init) for _ in 1:num_walkers]
    y_walkers = [copy(y_init) for _ in 1:num_walkers]
    weights   = ones(Float64, num_walkers)
    E_ref     = E_ref_initial
    E_history = Float64[]
    E_plot    = Float64[]

    # Pre-allocate
    drift_x_old = Vector{Vector{Float64}}(undef, num_walkers)
    drift_y_old = Vector{Vector{Float64}}(undef, num_walkers)
    drift_x_new = Vector{Vector{Float64}}(undef, num_walkers)
    drift_y_new = Vector{Vector{Float64}}(undef, num_walkers)
    E_loc_old   = Vector{Float64}(undef, num_walkers)
    E_loc_new   = Vector{Float64}(undef, num_walkers)
    new_x       = Vector{Vector{Float64}}(undef, num_walkers)
    new_y       = Vector{Vector{Float64}}(undef, num_walkers)

    # Compute initial drift forces and energies in one pass
    @threads for i in 1:num_walkers
        drift_x_old[i], drift_y_old[i], E_loc_old[i], _, _, _, _ = energy_estimators(
            x_walkers[i], y_walkers[i], L, R_match, Constants)
    end

    mode = quadratic ? "Quadratic (QDMC)" : "Linear DMC"
    prog = Progress(num_steps; desc="Running $mode...", showspeed=true)


    for step in 1:num_steps
        next!(prog)
        n = length(x_walkers)

        if length(new_x) != n
            resize!(new_x, n)
            resize!(new_y, n)
            resize!(E_loc_new, n)
            resize!(drift_x_new, n)
            resize!(drift_y_new, n)
        end

        if quadratic

            @threads for i in 1:n
                # Drift Δτ/2 at old position using precomputed drift
                x_d = wrap_position.(x_walkers[i] .+ drift_x_old[i] .* (Δτ/2), L)
                y_d = wrap_position.(y_walkers[i] .+ drift_y_old[i] .* (Δτ/2), L)

                # Diffusion
                x_d, y_d = diffusion_step(x_d, y_d, L, D, Δτ)

                # energy_estimators — single O(N^2) pass
                drift_x_new[i], drift_y_new[i], E_loc_new[i], _, _, _, _ = energy_estimators(
                    x_d, y_d, L, R_match, Constants)

                # Drift Δτ/2 at new position using new drift forces
                new_x[i] = wrap_position.(x_d .+ drift_x_new[i] .* (Δτ/2), L)
                new_y[i] = wrap_position.(y_d .+ drift_y_new[i] .* (Δτ/2), L)
            end

            x_walkers   = new_x[1:n]
            y_walkers   = new_y[1:n]
            drift_x_old = drift_x_new[1:n]
            drift_y_old = drift_y_new[1:n]

            avg_E = mean(E_loc_new[1:n])

            # # 5. Weights update
            # for i in 1:n
            #     weights[i] *= weight_update(E_loc_old[i], E_loc_new[i], E_ref, Δτ)
            # end

        else
            @threads for i in 1:n
                # Diffusion
                x_d, y_d = diffusion_step(x_walkers[i], y_walkers[i], L, D, Δτ)

                x_d = wrap_position.(x_d .+ drift_x_old[i] .* Δτ, L)
                y_d = wrap_position.(y_d .+ drift_y_old[i] .* Δτ, L)

                drift_x_old[i], drift_y_old[i], E_loc_new[i], _, _, _, _ = energy_estimators( x_d, y_d, L, R_match, Constants)                                       

                new_x[i] = x_d
                new_y[i] = y_d

            end

            x_walkers   = new_x[1:n]
            y_walkers   = new_y[1:n]

            avg_E = mean(E_loc_new[1:n])

            # # 5. Weights update
            # for i in 1:n
            #     weights[i] *= weight_update(E_loc_old[i], E_loc_new[i], E_ref, Δτ)
            # end
        end
        

        # # 6. Branching
        # num_copies = branching_step(weights[1:n])

        # new_x_b     = Vector{Vector{Float64}}()
        # new_y_b     = Vector{Vector{Float64}}()
        # new_dx      = Vector{Vector{Float64}}()
        # new_dy      = Vector{Vector{Float64}}()
        # E_loc_old_b = Float64[]

        # for i in 1:n
        
        #     for _ in 1:num_copies[i]
        #         push!(new_x_b, copy(x_walkers[i]))
        #         push!(new_y_b, copy(y_walkers[i]))
        #         push!(new_dx,  copy(drift_x_old[i]))
        #         push!(new_dy,  copy(drift_y_old[i]))
        #         push!(E_loc_old_b, E_loc_new[i])
        #     end
        # end

        # x_walkers   = new_x_b
        # y_walkers   = new_y_b
        # drift_x_old = new_dx
        # drift_y_old = new_dy
        # E_loc_old   = E_loc_old_b
        # weights     = ones(Float64, length(x_walkers))

        if isempty(x_walkers)
            @warn "All walkers died at step $step"
            break
        end

        # 7. Population control, with a 10% tolerance to prevent over-correction
        if abs(length(x_walkers)- num_target)/ num_target > 0.10
            E_ref = population_control(length(x_walkers), num_target, avg_E, Δτ)
        end

        # 8. Accumulate
        if step > num_equil
            push!(E_history, avg_E)
        end
        push!(E_plot, avg_E)
           if step % 10000 == 0
            println("\nStep $step: Walkers=$(length(x_walkers)), avg_E=$(round(avg_E/num_part, digits=5)), E_ref=$(round(E_ref/num_part, digits=5))")
        end
    end

    if plot_energy
        p_energy = plot(1:length(E_plot), E_plot ./ num_part,
             xlabel="DMC step", ylabel="E/N",
             title="DMC Energy History", linewidth=2)
        display(p_energy)
        println("\n>>> Press ENTER to continue...")
        readline()
    end

    isempty(E_history) && error("DMC failed: all walkers died before equilibration.")

    E_dmc     = mean(E_history)
    E_dmc_err = std(E_history) / sqrt(length(E_history))
    return E_dmc, E_dmc_err, E_history
end