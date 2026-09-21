# Stage1_2D_Dipol_System/scripts/dtau_convergence.jl
#
# Expects from main_DMC.jl: num_part, nr0_sq, L, R_opt, stage_dir,
# Δτ_DMC (from tune_delta_dmc.jl), num_walkers_dtau_study, τ_total,
# quadratic (Bool — which DMC scheme is being tested)

using Plots, LaTeXStrings

τ_factors = [0.05, 0.1, 0.5, 1.0, 2.0, 5.0]

vmc_run = load_run(result_path(stage_dir, "VMC", (N=num_part, nr0sq=nr0_sq)); run=1)
E_ref_initial = vmc_run.result.avg_energy

Constants = calculate_constants(L, R_opt)
trial = SingleLayerTrial(; R_match=R_opt, Constants=Constants)
coords_init = init_random_config((A=num_part,), L)

Δτ_vals, E_vals, E_err_vals = Float64[], Float64[], Float64[]
for factor in τ_factors
    Δτ = Δτ_DMC * factor
    num_steps_i = round(Int, τ_total / Δτ)
    num_equil_i = num_steps_i ÷ 5

    println("\n--- Δτ = $Δτ (factor $factor), num_steps = τ/Δτ = $num_steps_i ---")
    result = dmc(trial, coords_init, num_walkers_dtau_study, num_steps_i, Δτ, L,
                 E_ref_initial, num_walkers_dtau_study; num_equil=num_equil_i, quadratic=quadratic)
    push!(Δτ_vals, Δτ)
    push!(E_vals, result.E_dmc)
    push!(E_err_vals, result.E_dmc_err)
    println("E = $(result.E_dmc) ± $(result.E_dmc_err)")

    if plot_dmc_trace_diagnostic
        p_trace = plot_dmc_trace(result, num_equil_i; title="Δτ=$(round(Δτ,digits=6)), N=$(num_part), nr0²=$(nr0_sq)")
        display(p_trace)
        mkpath(joinpath(stage_dir, "data", "plots"))
        savefig(p_trace, joinpath(stage_dir, "data", "plots", "dmc_trace_N$(num_part)_nr0sq$(nr0_sq)_dtau$(round(Δτ,digits=6)).pdf"))
    end

    path = result_path(stage_dir, "dtau_convergence", (N=num_part, nr0sq=nr0_sq, factor=factor))
    save_run(path, (; num_part, nr0_sq, L, R_opt, Δτ, factor, num_steps_i, num_equil_i, τ_total,
                     num_walkers_dtau_study, quadratic), result; overwrite=force_rerun_overwrite)
end

intercept, intercept_err, coeffs, fitted_fn = weighted_extrapolation(
    Δτ_vals, E_vals, E_err_vals; powers = quadratic ? [2] : [1])

n_sigma = 1.0
Δτ_chosen = choose_largest_consistent(Δτ_vals, E_vals, E_err_vals, intercept, intercept_err; n_sigma=n_sigma)
println("\nΔτ_chosen = $Δτ_chosen  (largest Δτ within $(n_sigma)σ of the extrapolated Δτ→0 value)")

dtau_range = collect(LinRange(0, maximum(Δτ_vals), 100))
p = scatter(Δτ_vals, E_vals; yerror=E_err_vals, xlabel=L"\Delta\tau", ylabel=L"E",
            title="Δτ convergence, N=$(num_part), nr0²=$(nr0_sq)", label="DMC data", framestyle=:box)
plot!(dtau_range, fitted_fn.(dtau_range); label="fit ($(quadratic ? "quadratic" : "linear"))", linestyle=:dash)
scatter!([0.0], [intercept]; yerror=[intercept_err], label="extrapolated (Δτ→0)", marker=:star5, markersize=10, color=:red)
vline!([Δτ_chosen]; label="Δτ_chosen", linestyle=:dot, color=:black)
display(p)
mkpath(joinpath(stage_dir, "data", "plots"))
savefig(p, joinpath(stage_dir, "data", "plots", "dtau_convergence_N$(num_part)_nr0sq$(nr0_sq).pdf"))