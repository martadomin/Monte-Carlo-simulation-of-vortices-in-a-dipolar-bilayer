function move_one_part(coords::NamedTuple, delta::Float64, L::Float64)
    species, id = random_particle(coords)
    x_old, y_old = coords[species].x[id], coords[species].y[id]
    x_new = wrap_position(x_old + rand()*2*delta - delta, L)
    y_new = wrap_position(y_old + rand()*2*delta - delta, L)
    return species, id, x_new, y_new
end

function move_all_part(coords::NamedTuple, delta::Float64, L::Float64)
    return map(coords) do c
        (x = wrap_position.(c.x .+ randn(length(c.x)) .* delta, L),
         y = wrap_position.(c.y .+ randn(length(c.y)) .* delta, L))
    end
end

"""
    tune_delta(trial, coords, L; kwargs...) -> (delta, coords)

Tunes the Metropolis step size to achieve target acceptance ratio (~50%),
using single-particle moves and `Δlogpsi(trial, ...)`. Does not mutate
`coords` — works on and returns its own copy.
"""
function tune_delta(trial::TrialWavefunction, coords::NamedTuple, L::Float64;
    target_ratio::Float64 = 0.5, num_tune_steps::Int = 10^4, block_size::Int = 100)

    coords = copy_coords(coords)
    delta = 0.1 * L

    for _ in 1:num_tune_steps÷block_size
        accepted = 0
        for _ in 1:block_size
            species, id, x_new, y_new = move_one_part(coords, delta, L)
            Δlog_weight = 2 * Δlogpsi(trial, coords, species, id, x_new, y_new, L)
            if log(rand()) < Δlog_weight
                coords[species].x[id] = x_new
                coords[species].y[id] = y_new
                accepted += 1
            end
        end
        ratio = accepted / block_size
        delta *= (ratio + 1e-2) / target_ratio
        delta  = clamp(delta, 1e-6, L/2)
    end
    return delta, coords
end

"""
    metropolis(trial, counts, num_steps, delta, L; kwargs...) -> NamedTuple

Runs the Metropolis VMC loop for any TrialWavefunction, any number of
species. `counts` gives the particle count per species, e.g. `(A=30, B=30)`
for a two-species system or `(A=60,)` for one — see trial_wavefunction.jl
for the interface `trial` must implement.

# Returns
NamedTuple with fields: energies, energies_drift, energies_laplacian,
E_tot, E_sq, E_drift, E_laplacian, E_kin, E_int, acceptance_ratio,
coords, observables (stage-defined — see init_observables/record_observables!).
"""
function metropolis(trial::TrialWavefunction, counts::NamedTuple, num_steps::Int,
    delta::Float64, L::Float64;
    coords_init::Union{NamedTuple, Nothing} = nothing,
    final_energy_plot::Bool = false, plot_every::Int = 100,
    progress::Bool = true, num_bins::Int = 100, move_all::Bool = false)

    num_part = sum(values(counts))
    acceptance_count   = 0
    energies           = Float64[]
    energies_drift     = Float64[]
    energies_laplacian = Float64[]
    energy_plot        = Float64[]
    step_plot          = Int[]
    E_tot = E_sq = E_kin = E_int = E_drift_tot = E_laplacian_tot = 0.0
    n_uncorr = 0
    obs = init_observables(trial, num_bins)

    coords = isnothing(coords_init) ? init_random_config(counts, L) : copy_coords(coords_init)

    est = energy_estimators(trial, coords, L)
    @assert !isnan(est.E_total_local) "Initial configuration produced NaN energy. Re-initialize."
    E_local, E_local_drift, E_local_laplacian, E_kinetic, E_potential = est.E_total_local, est.E_total_drift, est.E_total_laplacian, est.E_kinetic, est.E_interaction
    log_weight_current = move_all ? logpsi(trial, coords, L) : 0.0
    progress_bar = progress ? Progress(num_steps; desc="Running Metropolis N=$num_part...", showspeed=true) : nothing

    for i in 1:num_steps
        progress && next!(progress_bar)

        if move_all
            coords_new = move_all_part(coords, delta, L)
            log_weight_new = logpsi(trial, coords_new, L)
            Δlog_weight = 2 * (log_weight_new - log_weight_current)
        else
            species, id, x_new, y_new = move_one_part(coords, delta, L)
            Δlog_weight = 2 * Δlogpsi(trial, coords, species, id, x_new, y_new, L)
        end

        if log(rand()) < Δlog_weight
            if move_all
                coords = coords_new
                log_weight_current = log_weight_new
            else
                coords[species].x[id] = x_new
                coords[species].y[id] = y_new
            end
            acceptance_count += 1
            est = energy_estimators(trial, coords, L)
            E_local, E_local_drift, E_local_laplacian, E_kinetic, E_potential =
                est.E_total_local, est.E_total_drift, est.E_total_laplacian, est.E_kinetic, est.E_interaction
        end

        if isnan(E_local)
            @warn "NaN detected at step $i"
            continue
        end

        push!(energies, E_local); push!(energies_drift, E_local_drift); push!(energies_laplacian, E_local_laplacian)
        E_tot += E_local; E_sq += E_local^2; E_kin += E_kinetic; E_int += E_potential
        E_drift_tot += E_local_drift; E_laplacian_tot += E_local_laplacian
        n_uncorr += 1

        record_observables!(trial, obs, coords, L)

        if i % plot_every == 0
            push!(step_plot, i); push!(energy_plot, E_local_drift / num_part)
        end
    end

    println("Acceptance ratio: ", acceptance_count / num_steps)

    if final_energy_plot
        p1 = plot(step_plot, energy_plot; xlabel="Step", ylabel="E/N (drift estimator)",
                  title="Energy evolution", legend=false, lw=2)
        mkpath(joinpath(stage_dir, "data", "plots"))
        display(p1)
        println("\n>>> Press ENTER to continue..."); readline()
    end

    return (; energies, energies_drift, energies_laplacian,
              E_tot=E_tot/n_uncorr, E_sq=E_sq/n_uncorr,
              E_drift=E_drift_tot/n_uncorr, E_laplacian=E_laplacian_tot/n_uncorr,
              E_kin=E_kin/n_uncorr, E_int=E_int/n_uncorr,
              acceptance_ratio=acceptance_count/num_steps,
              coords, observables=obs)
end