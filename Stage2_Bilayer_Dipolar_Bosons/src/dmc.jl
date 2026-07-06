using Statistics, Base.Threads

function diffusion_step(x_coord::Vector{Float64}, y_coord::Vector{Float64},
                        L::Float64, D::Float64, Δτ::Float64)::Tuple{Vector{Float64}, Vector{Float64}}
    x_proposed = wrap_position.(x_coord .+ sqrt(2 * D * Δτ) .* randn(length(x_coord)), L)
    y_proposed = wrap_position.(y_coord .+ sqrt(2 * D * Δτ) .* randn(length(y_coord)), L)
    return x_proposed, y_proposed
end

function weight_update(E_loc_old::Float64, E_loc_new::Float64,
                       E_ref::Float64, Δτ::Float64)::Float64
    if isnan(E_loc_old) || isnan(E_loc_new)
        return 0.0
    end
    exponent = -Δτ * (0.5*(E_loc_old + E_loc_new) - E_ref)
    return exp(clamp(exponent, -50.0, log(10.0)))
end

function branching_step(weights::Vector{Float64})::Vector{Int}
    safe_weights = [isnan(w) ? 0.0 : clamp(w, 0.0, 3.0) for w in weights]
    return [floor(Int, w + rand()) for w in safe_weights]
end

function population_control(num_walkers::Int, num_target::Int,
                             avg_E_loc::Float64, Δτ::Float64)::Float64
    α = 1.0
    return avg_E_loc - (α / Δτ) * log(num_walkers / num_target)
end

"""
    dmc(xA_init, yA_init, xB_init, yB_init,
        num_walkers, num_part, num_steps, Δτ,
        L, h, R_match, Constants, R0, itp_up, itp_upp,
        E_ref_initial, num_target;
        num_equil, quadratic, plot_energy)

Diffusion Monte Carlo for the bilayer dipolar system.
Supports both linear DMC and quadratic (2nd-order) DMC via the `quadratic` flag.

# Key differences from Stage I:
- Four coordinate arrays (xA, yA, xB, yB) instead of two
- `energy_estimators` requires h, R0, itp_up, itp_upp
- Drift forces are split into layer A (indices 1:N/2) and layer B (N/2+1:N)
  from the single `drift_x, drift_y` vector returned by energy_estimators
"""
function dmc(xA_init::Vector{Float64}, yA_init::Vector{Float64},
             xB_init::Vector{Float64}, yB_init::Vector{Float64},
             num_walkers::Int, num_part::Int, num_steps::Int,
             Δτ::Float64, L::Float64, h::Float64,
             R_match::Float64,
             Constants::Tuple{Float64, Float64, Float64},
             R0::Float64, itp_up, itp_upp,
             E_ref_initial::Float64, num_target::Int;
             num_equil::Int  = num_steps ÷ 5,
             quadratic::Bool = false,
             plot_energy::Bool = false)

    D      = 0.5
    N_half = num_part ÷ 2

    # ── Initialize walkers ───────────────────────────────────────────
    xA_walkers = [copy(xA_init) for _ in 1:num_walkers]
    yA_walkers = [copy(yA_init) for _ in 1:num_walkers]
    xB_walkers = [copy(xB_init) for _ in 1:num_walkers]
    yB_walkers = [copy(yB_init) for _ in 1:num_walkers]

    weights   = ones(Float64, num_walkers)
    E_ref     = E_ref_initial
    E_history = Float64[]
    E_plot    = Float64[]

    # ── Pre-allocate drift and energy arrays ─────────────────────────
    drift_xA_old = Vector{Vector{Float64}}(undef, num_walkers)
    drift_yA_old = Vector{Vector{Float64}}(undef, num_walkers)
    drift_xB_old = Vector{Vector{Float64}}(undef, num_walkers)
    drift_yB_old = Vector{Vector{Float64}}(undef, num_walkers)
    drift_xA_new = Vector{Vector{Float64}}(undef, num_walkers)
    drift_yA_new = Vector{Vector{Float64}}(undef, num_walkers)
    drift_xB_new = Vector{Vector{Float64}}(undef, num_walkers)
    drift_yB_new = Vector{Vector{Float64}}(undef, num_walkers)
    E_loc_old    = Vector{Float64}(undef, num_walkers)
    E_loc_new    = Vector{Float64}(undef, num_walkers)
    new_xA       = Vector{Vector{Float64}}(undef, num_walkers)
    new_yA       = Vector{Vector{Float64}}(undef, num_walkers)
    new_xB       = Vector{Vector{Float64}}(undef, num_walkers)
    new_yB       = Vector{Vector{Float64}}(undef, num_walkers)

    # ── Initial drift forces and energies ────────────────────────────
    @threads for i in 1:num_walkers
        dx, dy, E_loc_old[i], _, _, _, _ = energy_estimators(
            xA_walkers[i], yA_walkers[i], xB_walkers[i], yB_walkers[i],
            L, h, R_match, Constants, R0, itp_up, itp_upp)
        drift_xA_old[i] = dx[1:N_half]
        drift_yA_old[i] = dy[1:N_half]
        drift_xB_old[i] = dx[N_half+1:end]
        drift_yB_old[i] = dy[N_half+1:end]
    end

    mode = quadratic ? "Quadratic DMC" : "Linear DMC"
    println(mode)
    flush(stdout)
    prog = Progress(num_steps; desc="Running $mode...", showspeed=true)

    # ── Main DMC loop ────────────────────────────────────────────────
    for step in 1:num_steps
        next!(prog)
        n = length(xA_walkers)

        # Resize working arrays if population changed after branching
        for arr in (new_xA, new_yA, new_xB, new_yB,
                    E_loc_new, drift_xA_new, drift_yA_new, drift_xB_new, drift_yB_new)
            length(arr) != n && resize!(arr, n)
        end

        if quadratic
            # ── Quadratic (2nd-order) DMC ────────────────────────────
            @threads for i in 1:n

                # First half-drift
                x1A = wrap_position.(xA_walkers[i] .+ drift_xA_old[i] .* (Δτ/2), L)
                y1A = wrap_position.(yA_walkers[i] .+ drift_yA_old[i] .* (Δτ/2), L)
                x1B = wrap_position.(xB_walkers[i] .+ drift_xB_old[i] .* (Δτ/2), L)
                y1B = wrap_position.(yB_walkers[i] .+ drift_yB_old[i] .* (Δτ/2), L)

                F1x, F1y, _, _, _, _, _ = energy_estimators(
                    x1A, y1A, x1B, y1B, L, h, R_match, Constants, R0, itp_up, itp_upp)
                F1xA = F1x[1:N_half]; F1yA = F1y[1:N_half]
                F1xB = F1x[N_half+1:end]; F1yB = F1y[N_half+1:end]

                xdA = wrap_position.(xA_walkers[i] .+ 0.5*(drift_xA_old[i] .+ F1xA) .* (Δτ/2), L)
                ydA = wrap_position.(yA_walkers[i] .+ 0.5*(drift_yA_old[i] .+ F1yA) .* (Δτ/2), L)
                xdB = wrap_position.(xB_walkers[i] .+ 0.5*(drift_xB_old[i] .+ F1xB) .* (Δτ/2), L)
                ydB = wrap_position.(yB_walkers[i] .+ 0.5*(drift_yB_old[i] .+ F1yB) .* (Δτ/2), L)

                # Diffusion
                xdA, ydA = diffusion_step(xdA, ydA, L, D, Δτ)
                xdB, ydB = diffusion_step(xdB, ydB, L, D, Δτ)

                # Second half-drift
                Fmx, Fmy, _, _, _, _, _ = energy_estimators(
                    xdA, ydA, xdB, ydB, L, h, R_match, Constants, R0, itp_up, itp_upp)
                FmxA = Fmx[1:N_half]; FmyA = Fmy[1:N_half]
                FmxB = Fmx[N_half+1:end]; FmyB = Fmy[N_half+1:end]

                x2A = wrap_position.(xdA .+ FmxA .* (Δτ/2), L)
                y2A = wrap_position.(ydA .+ FmyA .* (Δτ/2), L)
                x2B = wrap_position.(xdB .+ FmxB .* (Δτ/2), L)
                y2B = wrap_position.(ydB .+ FmyB .* (Δτ/2), L)

                F2x, F2y, _, _, _, _, _ = energy_estimators(
                    x2A, y2A, x2B, y2B, L, h, R_match, Constants, R0, itp_up, itp_upp)
                F2xA = F2x[1:N_half]; F2yA = F2y[1:N_half]
                F2xB = F2x[N_half+1:end]; F2yB = F2y[N_half+1:end]

                new_xA[i] = wrap_position.(xdA .+ 0.5*(FmxA .+ F2xA) .* (Δτ/2), L)
                new_yA[i] = wrap_position.(ydA .+ 0.5*(FmyA .+ F2yA) .* (Δτ/2), L)
                new_xB[i] = wrap_position.(xdB .+ 0.5*(FmxB .+ F2xB) .* (Δτ/2), L)
                new_yB[i] = wrap_position.(ydB .+ 0.5*(FmyB .+ F2yB) .* (Δτ/2), L)

                dx, dy, E_loc_new[i], _, _, _, _ = energy_estimators(
                    new_xA[i], new_yA[i], new_xB[i], new_yB[i],
                    L, h, R_match, Constants, R0, itp_up, itp_upp)
                drift_xA_new[i] = dx[1:N_half]
                drift_yA_new[i] = dy[1:N_half]
                drift_xB_new[i] = dx[N_half+1:end]
                drift_yB_new[i] = dy[N_half+1:end]
            end

        else
            # ── Linear DMC ───────────────────────────────────────────
            @threads for i in 1:n

                # Diffusion + drift
                xdA, ydA = diffusion_step(xA_walkers[i], yA_walkers[i], L, D, Δτ)
                xdB, ydB = diffusion_step(xB_walkers[i], yB_walkers[i], L, D, Δτ)

                new_xA[i] = wrap_position.(xdA .+ drift_xA_old[i] .* Δτ, L)
                new_yA[i] = wrap_position.(ydA .+ drift_yA_old[i] .* Δτ, L)
                new_xB[i] = wrap_position.(xdB .+ drift_xB_old[i] .* Δτ, L)
                new_yB[i] = wrap_position.(ydB .+ drift_yB_old[i] .* Δτ, L)

                dx, dy, E_loc_new[i], _, _, _, _ = energy_estimators(
                    new_xA[i], new_yA[i], new_xB[i], new_yB[i],
                    L, h, R_match, Constants, R0, itp_up, itp_upp)
                drift_xA_new[i] = dx[1:N_half]
                drift_yA_new[i] = dy[1:N_half]
                drift_xB_new[i] = dx[N_half+1:end]
                drift_yB_new[i] = dy[N_half+1:end]
            end
        end

        # Update walker positions and drift forces
        xA_walkers   = new_xA[1:n]
        yA_walkers   = new_yA[1:n]
        xB_walkers   = new_xB[1:n]
        yB_walkers   = new_yB[1:n]
        drift_xA_old = drift_xA_new[1:n]
        drift_yA_old = drift_yA_new[1:n]
        drift_xB_old = drift_xB_new[1:n]
        drift_yB_old = drift_yB_new[1:n]

        avg_E = mean(E_loc_new[1:n])
        isnan(avg_E) && @warn "avg_E is NaN at step $step"

        # ── Weight update ────────────────────────────────────────────
        for i in 1:n
            weights[i] *= weight_update(E_loc_old[i], E_loc_new[i], E_ref, Δτ)
        end
        E_loc_old = copy(E_loc_new[1:n])

        # ── Branching ────────────────────────────────────────────────
        num_copies = branching_step(weights[1:n])

        new_xA_b   = Vector{Vector{Float64}}()
        new_yA_b   = Vector{Vector{Float64}}()
        new_xB_b   = Vector{Vector{Float64}}()
        new_yB_b   = Vector{Vector{Float64}}()
        new_dxA    = Vector{Vector{Float64}}()
        new_dyA    = Vector{Vector{Float64}}()
        new_dxB    = Vector{Vector{Float64}}()
        new_dyB    = Vector{Vector{Float64}}()
        E_loc_old_b = Float64[]

        for i in 1:n
            for _ in 1:num_copies[i]
                push!(new_xA_b,  copy(xA_walkers[i]))
                push!(new_yA_b,  copy(yA_walkers[i]))
                push!(new_xB_b,  copy(xB_walkers[i]))
                push!(new_yB_b,  copy(yB_walkers[i]))
                push!(new_dxA,   copy(drift_xA_old[i]))
                push!(new_dyA,   copy(drift_yA_old[i]))
                push!(new_dxB,   copy(drift_xB_old[i]))
                push!(new_dyB,   copy(drift_yB_old[i]))
                push!(E_loc_old_b, E_loc_new[i])
            end
        end

        isempty(new_xA_b) && (@warn "All walkers died at step $step"; break)

        xA_walkers   = new_xA_b
        yA_walkers   = new_yA_b
        xB_walkers   = new_xB_b
        yB_walkers   = new_yB_b
        drift_xA_old = new_dxA
        drift_yA_old = new_dyA
        drift_xB_old = new_dxB
        drift_yB_old = new_dyB
        E_loc_old    = E_loc_old_b
        weights      = ones(Float64, length(xA_walkers))

        # ── Population control ───────────────────────────────────────
        avg_E_post = isempty(E_loc_old_b) ? E_ref_initial : mean(E_loc_old_b)
        E_ref      = population_control(length(xA_walkers), num_target, avg_E_post, Δτ)

        # ── Accumulate after equilibration ───────────────────────────
        step > num_equil && push!(E_history, avg_E_post)
        push!(E_plot, avg_E)

        if step % 10000 == 0
            println("\nStep $step: walkers=$(length(xA_walkers)), " *
                    "E/N=$(round(avg_E/num_part, digits=5)), " *
                    "E_ref/N=$(round(E_ref/num_part, digits=5))")
            flush(stdout)
        end
    end

    if plot_energy
        p_e = plot(1:length(E_plot), E_plot ./ num_part;
                   xlabel="DMC step", ylabel="E/N",
                   title="DMC Energy History", linewidth=2, legend=false)
        display(p_e)
        println("\n>>> Press ENTER to continue...")
        readline()
    end

    isempty(E_history) && error("DMC failed: all walkers died before equilibration.")

    E_dmc     = mean(E_history)
    E_dmc_err = std(E_history) / sqrt(length(E_history))
    return E_dmc, E_dmc_err, E_history
end