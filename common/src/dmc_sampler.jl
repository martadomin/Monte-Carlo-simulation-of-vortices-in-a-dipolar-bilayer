# common/src/dmc.jl
#
# Generic DMC driver for any TrialWavefunction/species count. Dispatches
# on `trial` for energy_estimators, init_observables, record_observables!
# (see trial_wavefunction.jl). Uses dmc_kernel.jl's diffusion_step,
# weight_update, branching_step, population_control — none of which know
# about species, only about walker energies/weights as plain numbers.

using Statistics, Base.Threads, ProgressMeter, Plots

"""
    dmc(trial, coords_init, num_walkers, num_steps, Δτ, L,
        E_ref_initial, num_target; kwargs...) -> NamedTuple

Diffusion Monte Carlo for any TrialWavefunction. Supports linear and
quadratic (2nd-order) DMC via `quadratic`.

# Density / g(r) accumulation (mixed estimator)
Observables are accumulated from the POST-BRANCHING walker population,
once per step past equilibration: after `branching_step`, a walker that
produced k copies literally appears k times in `walkers`, so looping once
over every entry already carries the correct branching weight.

# Returns
NamedTuple with fields: E_dmc, E_dmc_err, E_history, observables (stage-
defined — see init_observables/record_observables!).
"""
function dmc(trial::TrialWavefunction, coords_init::NamedTuple,
             num_walkers::Int, num_steps::Int, Δτ::Float64, L::Float64,
             E_ref_initial::Float64, num_target::Int;
             num_equil::Int = num_steps ÷ 5, quadratic::Bool = false,
             plot_energy::Bool = false, num_bins::Int = 100)

    D = 0.5
    walkers   = [copy_coords(coords_init) for _ in 1:num_walkers]
    weights   = ones(Float64, num_walkers)
    E_ref     = E_ref_initial
    E_history = Float64[]
    E_plot    = Float64[]
    obs       = init_observables(trial, num_bins)
    n_obs_samples = 0

    drift_old = Vector{NamedTuple}(undef, num_walkers)
    E_loc_old = Vector{Float64}(undef, num_walkers)

    @threads for i in 1:num_walkers
        est = energy_estimators(trial, walkers[i], L)
        drift_old[i] = est.drift
        E_loc_old[i] = est.E_total_local
    end

    mode = quadratic ? "Quadratic DMC" : "Linear DMC"
    println(mode); flush(stdout)
    prog = Progress(num_steps; desc="Running $mode...", showspeed=true)

    for step in 1:num_steps
        next!(prog)
        n = length(walkers)
        new_walkers = Vector{NamedTuple}(undef, n)
        drift_new   = Vector{NamedTuple}(undef, n)
        E_loc_new   = Vector{Float64}(undef, n)

        @threads for i in 1:n
            if quadratic
                half1 = map(walkers[i], drift_old[i]) do c, d
                    (x = wrap_position.(c.x .+ d.x .* (Δτ/2), L),
                     y = wrap_position.(c.y .+ d.y .* (Δτ/2), L))
                end
                est1 = energy_estimators(trial, half1, L)

                mid = map(walkers[i], drift_old[i], est1.drift) do c, d0, d1
                    (x = wrap_position.(c.x .+ 0.5*(d0.x .+ d1.x) .* (Δτ/2), L),
                     y = wrap_position.(c.y .+ 0.5*(d0.y .+ d1.y) .* (Δτ/2), L))
                end
                mid = diffusion_step(mid, L, D, Δτ)
                estm = energy_estimators(trial, mid, L)

                half2 = map(mid, estm.drift) do c, d
                    (x = wrap_position.(c.x .+ d.x .* (Δτ/2), L),
                     y = wrap_position.(c.y .+ d.y .* (Δτ/2), L))
                end
                est2 = energy_estimators(trial, half2, L)

                new_walkers[i] = map(mid, estm.drift, est2.drift) do c, dm, d2
                    (x = wrap_position.(c.x .+ 0.5*(dm.x .+ d2.x) .* (Δτ/2), L),
                     y = wrap_position.(c.y .+ 0.5*(dm.y .+ d2.y) .* (Δτ/2), L))
                end
                est_final = energy_estimators(trial, new_walkers[i], L)
            else
                diffused = diffusion_step(walkers[i], L, D, Δτ)
                new_walkers[i] = map(diffused, drift_old[i]) do c, d
                    (x = wrap_position.(c.x .+ d.x .* Δτ, L),
                     y = wrap_position.(c.y .+ d.y .* Δτ, L))
                end
                est_final = energy_estimators(trial, new_walkers[i], L)
            end
            drift_new[i] = est_final.drift
            E_loc_new[i] = est_final.E_total_local
        end

        walkers, drift_old = new_walkers, drift_new
        avg_E = mean(E_loc_new)
        isnan(avg_E) && @warn "avg_E is NaN at step $step"

        for i in 1:n
            weights[i] *= weight_update(E_loc_old[i], E_loc_new[i], E_ref, Δτ)
        end
        E_loc_old = copy(E_loc_new)

        num_copies = branching_step(weights[1:n])
        branched_walkers, branched_drift, branched_E = NamedTuple[], NamedTuple[], Float64[]
        for i in 1:n, _ in 1:num_copies[i]
            push!(branched_walkers, copy_coords(walkers[i]))
            push!(branched_drift, map(d -> (x = copy(d.x), y = copy(d.y)), drift_old[i]))
            push!(branched_E, E_loc_new[i])
        end

        isempty(branched_walkers) && (@warn "All walkers died at step $step"; break)

        walkers, drift_old, E_loc_old = branched_walkers, branched_drift, branched_E
        weights = ones(Float64, length(walkers))

        avg_E_post = mean(E_loc_old)
        E_ref = population_control(length(walkers), num_target, avg_E_post, Δτ)

        if step > num_equil
            for w in walkers
                record_observables!(trial, obs, w, L)
            end
            n_obs_samples += length(walkers)
            push!(E_history, avg_E_post)
        end
        push!(E_plot, avg_E)

        if step % 10000 == 0
            println("\nStep $step: walkers=$(length(walkers)), " *
                    "E/N=$(round(avg_E/num_particles(walkers[1]), digits=5)), " *
                    "E_ref/N=$(round(E_ref/num_particles(walkers[1]), digits=5))")
            flush(stdout)
        end
    end

    if plot_energy
        p_e = plot(1:length(E_plot), E_plot ./ num_particles(walkers[1]);
                   xlabel="DMC step", ylabel="E/N", title="DMC Energy History",
                   linewidth=2, legend=false)
        display(p_e)
        println("\n>>> Press ENTER to continue..."); readline()
    end

    isempty(E_history) && error("DMC failed: all walkers died before equilibration.")
    @assert n_obs_samples > 0 "No observable samples accumulated (num_equil >= num_steps?)."

    return (; E_dmc = mean(E_history), E_dmc_err = std(E_history)/sqrt(length(E_history)),
            E_history, E_plot, observables = obs)
end