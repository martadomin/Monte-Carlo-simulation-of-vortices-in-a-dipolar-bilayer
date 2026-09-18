# Stage2_Bilayer_Dipolar_Bosons/scripts/dtau_convergence.jl
#
# Expects from main_DMC.jl: num_part, nr0_sq, h, L, r_min, Δ_shoot, tol_shoot,
# Δτ_DMC (from tune_delta_dmc.jl), num_walkers_dtau_study, num_steps_dmc,
# quadratic, stage_dir, stage1_dir

using Plots, LaTeXStrings

τ_factors = [0.05, 0.1, 0.5, 1.0, 2.0, 5.0]

N_half = num_part ÷ 2
R_match = load_Rmatch(stage1_dir, N_half, L)
Constants = calculate_constants(L, R_match)

R0_run = load_run(result_path(stage_dir, "R0_optimum", (N=num_part, nr0sq=nr0_sq, h=h)); run=1)
R0_opt = R0_run.result.R0_opt

fAB_result = build_fAB(h, R0_opt, r_min, Δ_shoot, tol_shoot, nr0_sq, num_part)
fAB_result === nothing && error("build_fAB failed at R0_opt = $R0_opt — cannot run Δτ convergence.")
_, _, itp_u, itp_up, itp_upp, _ = fAB_result

trial = BilayerTrial(; R_match, Constants, R0=R0_opt, itp_u, itp_up, itp_upp, h)

vmc_run = load_run(result_path(stage_dir, "VMC", (N=num_part, nr0sq=nr0_sq, h=h)); run=1)
E_ref_initial = vmc_run.result.avg_energy

coords_init = init_random_config((A=N_half, B=N_half), L)

Δτ_vals, E_vals, E_err_vals = Float64[], Float64[], Float64[]
for factor in τ_factors
    Δτ = Δτ_DMC * factor
    println("\n--- Δτ = $Δτ (factor $factor) ---")
    result = dmc(trial, coords_init, num_walkers_dtau_study, num_steps_dmc, Δτ, L,
                 E_ref_initial, num_walkers_dtau_study; num_equil=num_steps_dmc ÷ 5, quadratic=quadratic)
    push!(Δτ_vals, Δτ)
    push!(E_vals, result.E_dmc)
    push!(E_err_vals, result.E_dmc_err)
    println("E = $(result.E_dmc) ± $(result.E_dmc_err)")

    if plot_dmc_trace_diagnostic
        p_trace = plot_dmc_trace(result, num_steps_dmc ÷ 5; title="Δτ=$(round(Δτ,digits=6)), N=$(num_part), nr0²=$(nr0_sq), h=$(h)")
        display(p_trace)
        savefig(p_trace, joinpath(stage_dir, "data", "plots", "dmc_trace_N$(num_part)_nr0sq$(nr0_sq)_h$(h)_dtau$(round(Δτ,digits=6)).pdf"))
    end

    path = result_path(stage_dir, "dtau_convergence", (N=num_part, nr0sq=nr0_sq, h=h, factor=factor))
    save_run(path, (; num_part, nr0_sq, h, L, R_match, R0=R0_opt, Δτ, factor,
                     num_walkers_dtau_study, num_steps_dmc, quadratic), result)
end

intercept, intercept_err, coeffs, fitted_fn = weighted_extrapolation(
    Δτ_vals, E_vals, E_err_vals; powers = quadratic ? [2] : [1])

n_sigma = 1.0
Δτ_chosen = choose_largest_consistent(Δτ_vals, E_vals, E_err_vals, intercept, intercept_err; n_sigma=n_sigma)
println("\nΔτ_chosen = $Δτ_chosen  (largest Δτ within $(n_sigma)σ of the extrapolated Δτ→0 value)")

dtau_range = collect(LinRange(0, maximum(Δτ_vals), 100))
p = scatter(Δτ_vals, E_vals; yerror=E_err_vals, xlabel=L"\Delta\tau", ylabel=L"E",
            title="Δτ convergence, N=$(num_part), nr0²=$(nr0_sq), h=$(h)", label="DMC data", framestyle=:box)
plot!(dtau_range, fitted_fn.(dtau_range); label="fit ($(quadratic ? "quadratic" : "linear"))", linestyle=:dash)
scatter!([0.0], [intercept]; yerror=[intercept_err], label="extrapolated (Δτ→0)", marker=:star5, markersize=10, color=:red)
vline!([Δτ_chosen]; label="Δτ_chosen", linestyle=:dot, color=:black)
display(p)
savefig(p, joinpath(stage_dir, "data", "plots", "dtau_convergence_N$(num_part)_nr0sq$(nr0_sq)_h$(h).pdf"))