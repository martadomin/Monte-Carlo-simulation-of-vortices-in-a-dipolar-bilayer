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
        moved_id          = rand(1:length(x_A))
        x_A_new[moved_id] = wrap_position(x_A[moved_id] + rand() * 2 * delta - delta, L)
        y_A_new[moved_id] = wrap_position(y_A[moved_id] + rand() * 2 * delta - delta, L)
        return moved_id, :A, x_A_new, y_A_new, x_B_new, y_B_new
    else
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

"""
    compute_logΨ(x_A, y_A, x_B, y_B, R_match, R0, itp_u, L, Constants) -> Float64

Computes log|Ψ_T| for the bilayer trial wavefunction:
    log Ψ = Σ_{i<j} u_AA(r_ij) + Σ_{α<β} u_AA(r_αβ) + Σ_{i,α} u_AB(r_iα)
Used in the move_all Metropolis acceptance ratio.
"""
function compute_logΨ(x_A::Vector{Float64}, y_A::Vector{Float64},
                      x_B::Vector{Float64}, y_B::Vector{Float64},
                      R_match::Float64, R0::Float64, itp_u,
                      L::Float64, Constants::Tuple{Float64, Float64, Float64},
                      x_vortex_A::Float64, y_vortex_A::Float64,
                      x_vortex_B::Float64, y_vortex_B::Float64,
                      lA::Float64, lB::Float64)::Float64
    logΨ   = 0.0
    N_half = length(x_A)

    # AA pairs
    @inbounds for i in 1:N_half
        @inbounds for j in (i + 1):N_half
            dx = get_periodic_difference(x_A[i], x_A[j], L)
            dy = get_periodic_difference(y_A[i], y_A[j], L)
            r  = sqrt(dx^2 + dy^2)
            logΨ += u_AA(r, R_match, L, Constants)
        end
    end

    # BB pairs
    @inbounds for α in 1:N_half
        @inbounds for β in (α + 1):N_half
            dx = get_periodic_difference(x_B[α], x_B[β], L)
            dy = get_periodic_difference(y_B[α], y_B[β], L)
            r  = sqrt(dx^2 + dy^2)
            logΨ += u_AA(r, R_match, L, Constants)
        end
    end

    # AB pairs
    @inbounds for i in 1:N_half
        @inbounds for α in 1:N_half
            dx = get_periodic_difference(x_A[i], x_B[α], L)
            dy = get_periodic_difference(y_A[i], y_B[α], L)
            r  = sqrt(dx^2 + dy^2)
            logΨ += u_AB(r, R0, itp_u)
        end
    end

    # Vortex one-body terms — one per particle, NOT pairwise
    @inbounds for i in 1:N_half
        rA = sqrt(get_periodic_difference(x_A[i], x_vortex_A, L)^2 +
                  get_periodic_difference(y_A[i], y_vortex_A, L)^2)
        logΨ += u_vortex(rA, lA, L)
    end
    @inbounds for α in 1:N_half
        rB = sqrt(get_periodic_difference(x_B[α], x_vortex_B, L)^2 +
                  get_periodic_difference(y_B[α], y_vortex_B, L)^2)
        logΨ += u_vortex(rB, lB, L)
    end

    return logΨ
end
"""
    compute_ΔlogΨ(...) -> Float64

Computes the change in log|Ψ_T| when one particle is moved.
Only valid for single-particle moves (move_one_part).
"""
function compute_ΔlogΨ(x_A_old::Vector{Float64}, y_A_old::Vector{Float64},
                        x_B_old::Vector{Float64}, y_B_old::Vector{Float64},
                        x_A_new::Vector{Float64}, y_A_new::Vector{Float64},
                        x_B_new::Vector{Float64}, y_B_new::Vector{Float64},
                        id::Int, layer::Symbol,
                        R_match::Float64, L::Float64,
                        Constants::Tuple{Float64, Float64, Float64},
                        R0::Float64, itp_u,
                        x_vortex_A::Float64, y_vortex_A::Float64,
                        x_vortex_B::Float64, y_vortex_B::Float64,
                        lA::Float64, lB::Float64)::Float64

    ΔlogΨ  = 0.0
    N_half = length(x_A_old)

    if layer == :A
        @inbounds for j in 1:N_half
            if j != id
                r_old = sqrt(get_periodic_difference(x_A_old[id], x_A_old[j], L)^2 +
                             get_periodic_difference(y_A_old[id], y_A_old[j], L)^2)
                r_new = sqrt(get_periodic_difference(x_A_new[id], x_A_old[j], L)^2 +
                             get_periodic_difference(y_A_new[id], y_A_old[j], L)^2)
                ΔlogΨ += u_AA(r_new, R_match, L, Constants) -
                         u_AA(r_old, R_match, L, Constants)
            end
            r_old = sqrt(get_periodic_difference(x_A_old[id], x_B_old[j], L)^2 +
                         get_periodic_difference(y_A_old[id], y_B_old[j], L)^2)
            r_new = sqrt(get_periodic_difference(x_A_new[id], x_B_old[j], L)^2 +
                         get_periodic_difference(y_A_new[id], y_B_old[j], L)^2)
            ΔlogΨ += u_AB(r_new, R0, itp_u) - u_AB(r_old, R0, itp_u)
        end

        # Vortex one-body term for the moved A particle
        rA_old = sqrt(get_periodic_difference(x_A_old[id], x_vortex_A, L)^2 +
                      get_periodic_difference(y_A_old[id], y_vortex_A, L)^2)
        rA_new = sqrt(get_periodic_difference(x_A_new[id], x_vortex_A, L)^2 +
                      get_periodic_difference(y_A_new[id], y_vortex_A, L)^2)
        ΔlogΨ += u_vortex(rA_new, lA, L) - u_vortex(rA_old, lA, L)

    else  # layer == :B
        @inbounds for j in 1:N_half
            if j != id
                r_old = sqrt(get_periodic_difference(x_B_old[id], x_B_old[j], L)^2 +
                             get_periodic_difference(y_B_old[id], y_B_old[j], L)^2)
                r_new = sqrt(get_periodic_difference(x_B_new[id], x_B_old[j], L)^2 +
                             get_periodic_difference(y_B_new[id], y_B_old[j], L)^2)
                ΔlogΨ += u_AA(r_new, R_match, L, Constants) -
                         u_AA(r_old, R_match, L, Constants)
            end
            r_old = sqrt(get_periodic_difference(x_B_old[id], x_A_old[j], L)^2 +
                         get_periodic_difference(y_B_old[id], y_A_old[j], L)^2)
            r_new = sqrt(get_periodic_difference(x_B_new[id], x_A_old[j], L)^2 +
                         get_periodic_difference(y_B_new[id], y_A_old[j], L)^2)
            ΔlogΨ += u_AB(r_new, R0, itp_u) - u_AB(r_old, R0, itp_u)
        end

        # Vortex one-body term for the moved B particle
        rB_old = sqrt(get_periodic_difference(x_B_old[id], x_vortex_B, L)^2 +
                      get_periodic_difference(y_B_old[id], y_vortex_B, L)^2)
        rB_new = sqrt(get_periodic_difference(x_B_new[id], x_vortex_B, L)^2 +
                      get_periodic_difference(y_B_new[id], y_vortex_B, L)^2)
        ΔlogΨ += u_vortex(rB_new, lB, L) - u_vortex(rB_old, lB, L)
    end

    return ΔlogΨ
end
"""
    tune_delta(...) -> (delta, x_A, y_A, x_B, y_B)

Tunes the Metropolis step size to achieve target acceptance ratio (~50%).
Uses single-particle moves only.
"""
function tune_delta(
    x_A::Vector{Float64},
    y_A::Vector{Float64},
    x_B::Vector{Float64},
    y_B::Vector{Float64},
    L::Float64,
    R_match::Float64,
    Constants::Tuple{Float64, Float64, Float64},
    x_vortex_A::Float64, y_vortex_A::Float64,
    x_vortex_B::Float64, y_vortex_B::Float64,
    lA::Float64, lB::Float64,
    R0::Float64,
    itp_u;
    target_ratio::Float64 = 0.5,
    num_tune_steps::Int   = 10^4,
    block_size::Int       = 100
)::Tuple{Float64, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}}

    delta = 0.1 * L

    for _ in 1:num_tune_steps÷block_size
        accepted = 0
        for _ in 1:block_size
            moved_id, layer, x_A_new, y_A_new, x_B_new, y_B_new =
            move_one_part(x_A, y_A, x_B, y_B, delta, L)
            ΔlogΨ = compute_ΔlogΨ(x_A, y_A, x_B, y_B,
                                x_A_new, y_A_new, x_B_new, y_B_new,
                                moved_id, layer, R_match, L, Constants, R0, itp_u,
                                x_vortex_A, y_vortex_A, x_vortex_B, y_vortex_B, lA, lB)
            if log(rand()) < 2 * ΔlogΨ
                x_A = x_A_new
                y_A = y_A_new
                x_B = x_B_new
                y_B = y_B_new
                accepted += 1
            end
        end
        ratio = accepted / block_size
        delta *= (ratio + 1e-2) / target_ratio
        delta  = min(delta, L/2)
        delta  = max(delta, 1e-6)
    end

    return delta, x_A, y_A, x_B, y_B
end

"""
    metropolis(num_part, num_steps, delta, L, h, R_match, R0, itp_u, itp_up, itp_upp,
               Constants; kwargs...) -> (energies, ..., x_coord, y_coord)

Runs the Metropolis VMC loop for the bilayer dipolar system.

# Input:
- `num_part`   : Total number of particles (N = N_A + N_B, must be even).
- `num_steps`  : Number of MC steps.
- `delta`      : Step size (tuned externally via tune_delta).
- `L`          : Box length.
- `h`          : Interlayer separation.
- `R_match`    : Matching radius for the intra-layer Jastrow factor.
- `R0`         : Cutoff radius for the inter-layer Jastrow factor (variational parameter).
- `itp_u`      : Cubic spline interpolant of log(f_AB), from build_fAB.
- `itp_up`     : Cubic spline interpolant of u'_AB, from build_fAB.
- `itp_upp`    : Interpolant of u''_AB, from build_fAB.
- `Constants`  : (C1, C2, C3) for the intra-layer Jastrow factor.

# Keyword arguments:
- `x_A_init, y_A_init, x_B_init, y_B_init` : Initial positions (optional).
- `final_energy_plot` : Show energy trace plot at the end.
- `plot_every`        : Stride for energy trace accumulation.
- `progress`          : Show progress bar.
- `move_all`          : Use move_all_part instead of move_one_part.
"""
function metropolis(
    num_part::Int,
    num_steps::Int,
    delta::Float64,
    L::Float64,
    h::Float64,
    R_match::Float64,
    R0::Float64,
    itp_u,
    itp_up,
    itp_upp,
    Constants::Tuple{Float64, Float64, Float64};
    x_vortex_A::Float64, y_vortex_A::Float64,
    x_vortex_B::Float64, y_vortex_B::Float64,
    lA::Float64,
    lB::Float64,
    x_A_init::Union{Vector{Float64}, Nothing} = nothing,
    y_A_init::Union{Vector{Float64}, Nothing} = nothing,
    x_B_init::Union{Vector{Float64}, Nothing} = nothing,
    y_B_init::Union{Vector{Float64}, Nothing} = nothing,
    final_energy_plot::Bool = false,
    plot_every::Int         = 100,
    progress::Bool          = true,
    num_bins::Int           = 100,
    move_all::Bool          = false
)::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64},
         Float64, Float64, Float64, Float64, Float64, Float64, Float64,
         Vector{Float64}, Vector{Float64},
         Matrix{Float64}, Matrix{Float64}, Vector{Float64},
         Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, 
         Vector{Float64}}
         
    N_half = num_part ÷ 2

    acceptance_ratio  = 0.0
    n_uncorr          = 0
    step_block        = 1
    E_tot             = 0.0
    E_sq              = 0.0
    E_kin             = 0.0
    E_int             = 0.0
    E_tot_drift       = 0.0
    E_tot_laplacian   = 0.0
    energies          = Float64[]
    energies_drift    = Float64[]
    energies_laplacian= Float64[]
    energy_plot       = Float64[]
    step_plot         = Int[]
    gAA_r               = zeros(Float64, num_bins)
    gBB_r               = zeros(Float64, num_bins)
    gAB_r               = zeros(Float64, num_bins)
    n_gr_samples      = 0
    E_local           = NaN
    E_local_drift     = NaN
    E_local_laplacian = NaN
    E_kinetic         = NaN
    E_potential       = NaN
    n_xy_A            = zeros(Float64, num_bins, num_bins)
    n_xy_B            = zeros(Float64, num_bins, num_bins)

    # ── Initial configuration ────────────────────────────────────────
    if any(isnothing, (x_A_init, y_A_init, x_B_init, y_B_init))
        x_coord, y_coord = random_initial_config(num_part, L, "Uniform")
        x_A = x_coord[1:N_half]
        y_A = y_coord[1:N_half]
        x_B = x_coord[N_half+1:end]
        y_B = y_coord[N_half+1:end]
    else
        x_A = copy(x_A_init)
        y_A = copy(y_A_init)
        x_B = copy(x_B_init)
        y_B = copy(y_B_init)
    end

    # ── Initial energy ───────────────────────────────────────────────
    _, _, E_local, E_local_drift, E_local_laplacian, E_kinetic, E_potential = energy_estimators(x_A, y_A, x_B, y_B,
                                                                                                x_vortex_A, y_vortex_A, x_vortex_B, y_vortex_B,
                                                                                                L, h, lA, lB, R_match, Constants, R0, itp_up, itp_upp)
    @assert !isnan(E_local) "Initial configuration produced NaN energy. Re-initialize."

    if progress
        progress_bar = Progress(num_steps; desc="Running Metropolis N=$num_part...", showspeed=true)
    end

    logΨ_current = move_all ? compute_logΨ(x_A, y_A, x_B, y_B, R_match, R0, itp_u, L, Constants,
                                            x_vortex_A, y_vortex_A, x_vortex_B, y_vortex_B, lA, lB) : 0.0

    # ── Main MC loop ─────────────────────────────────────────────────
    for i in 1:num_steps
        progress && next!(progress_bar)

        if move_all
            x_A_new, y_A_new, x_B_new, y_B_new =
                move_all_part(x_A, y_A, x_B, y_B, delta, L)
            logΨ_new = compute_logΨ(x_A_new, y_A_new, x_B_new, y_B_new,
                                    R_match, R0, itp_u, L, Constants,
                                    x_vortex_A, y_vortex_A, x_vortex_B, y_vortex_B, lA, lB)
            ΔlogΨ = 2 * (logΨ_new - logΨ_current)
        else
        moved_id, layer, x_A_new, y_A_new, x_B_new, y_B_new = move_one_part(x_A, y_A, x_B, y_B, delta, L)
        ΔlogΨ = 2 * compute_ΔlogΨ(x_A, y_A, x_B, y_B,
                                x_A_new, y_A_new, x_B_new, y_B_new,
                                moved_id, layer,
                                R_match, L, Constants, R0, itp_u,
                                x_vortex_A, y_vortex_A, x_vortex_B, y_vortex_B, lA, lB)
        end

        if log(rand()) < ΔlogΨ
            x_A = x_A_new
            y_A = y_A_new
            x_B = x_B_new
            y_B = y_B_new
            acceptance_ratio += 1
            move_all && (logΨ_current = logΨ_new)

            _, _, E_local, E_local_drift, E_local_laplacian, E_kinetic, E_potential =
                energy_estimators(x_A, y_A, x_B, y_B, x_vortex_A, y_vortex_A, x_vortex_B, y_vortex_B, L, h, lA, lB, R_match, Constants, R0, itp_up, itp_upp)
        end

        if i % step_block == 0
            if isnan(E_local)
                @warn "NaN detected at step $i"
                continue
            end

            push!(energies,           E_local)
            push!(energies_drift,     E_local_drift)
            push!(energies_laplacian, E_local_laplacian)

            E_tot           += E_local
            E_sq            += E_local^2
            E_kin           += E_kinetic
            E_int           += E_potential
            E_tot_drift     += E_local_drift
            E_tot_laplacian += E_local_laplacian
            n_uncorr        += 1

            # Accumulate n_xy_A and n_xy_B
            accumulate_density!(n_xy_A, x_A, y_A, L)
            accumulate_density!(n_xy_B, x_B, y_B, L)

            # Accumulate gr
            accumulate_gr!(gAA_r, x_A, y_A, L)
            accumulate_gr!(gBB_r, x_B, y_B, L)
            accumulate_g_AB_r!(gAB_r, x_A, y_A, x_B, y_B, L)

            if i % plot_every == 0
                push!(step_plot,   i)
                push!(energy_plot, E_local_drift / num_part)
            end
            n_gr_samples += 1
        end
    end

    println("Acceptance ratio: ", acceptance_ratio / num_steps)

    if final_energy_plot
        p1 = plot(step_plot, energy_plot;
                  xlabel  = "Step",
                  ylabel  = "E/N (drift estimator)",
                  title   = "Energy evolution",
                  legend  = false,
                  lw      = 2)
        display(p1)
        println("\n>>> Press ENTER to continue...")
        readline()
    end

    x_coord = vcat(x_A, x_B)
    y_coord = vcat(y_A, y_B)

    g_total_r = gAA_r .+ gBB_r .+ 2 .* gAB_r

    xy_bins = normalize_density!(n_xy_A, n_uncorr, L)
    normalize_density!(n_xy_B, n_uncorr, L)

    r_vals = normalize_gr!(gAA_r, N_half, L, n_gr_samples)
    normalize_gr!(gBB_r, N_half, L, n_gr_samples)
    normalize_gr!(gAB_r, N_half, L, n_gr_samples)
    normalize_gr!(g_total_r, num_part, L, n_gr_samples)

    return energies, energies_drift, energies_laplacian,
           E_tot/n_uncorr, E_sq/n_uncorr, E_tot_drift/n_uncorr, E_tot_laplacian/n_uncorr,
           E_kin/n_uncorr, E_int/n_uncorr, acceptance_ratio/num_steps,
           x_coord, y_coord,
           n_xy_A, n_xy_B, xy_bins,
           gAA_r, gBB_r, gAB_r, g_total_r, r_vals
end