# Stage1_2D_Dipol_System/scripts/num_walkers_convergence.jl
#
# Expects from main_DMC.jl: num_part, nr0_sq, L, R_opt, stage_dir,
# Δτ_chosen (from dtau_convergence.jl), num_walkers_vals, τ_total,
# quadratic

using Plots, LaTeXStrings

num_steps = round(Int, τ_total / Δτ_chosen)
num_equil = num_steps ÷ 5

println("Using num_steps = τ_total/Δτ_chosen = $num_steps (same τ_total=$τ_total as the Δτ study)")

Constants = calculate_constants(L, R_opt)
trial = SingleLayerTrial(; R_match=R_opt, Constants=Constants)

vmc_run = load_run(result_path(stage_dir, "VMC", (N=num_part, nr0sq=nr0_sq)); run=1)
E_ref_initial = vmc_run.result.avg_energy

coords_init = init_random_config((A=num_part,), L)

E_vals, E_err_vals = Float64[], Float64[]
for num_walkers in num_walkers_vals
    println("\n--- num_walkers = $num_walkers ---")
    result = dmc(trial, coords_init, num_walkers, num_steps, Δτ_chosen, L,
                 E_ref_initial, num_walkers; num_equil=num_equil, quadratic=quadratic)
    push!(E_vals, result.E_dmc)
    push!(E_err_vals, result.E_dmc_err)
    println("E = $(result.E_dmc) ± $(result.E_dmc_err)")

    if plot_dmc_trace_diagnostic
        p_trace = plot_dmc_trace(result, num_equil; title="num_walkers=$(num_walkers), N=$(num_part), nr0²=$(nr0_sq)")
        display(p_trace)
        mkpath(joinpath(stage_dir, "data", "plots"))
        savefig(p_trace, joinpath(stage_dir, "data", "plots", "dmc_trace_N$(num_part)_nr0sq$(nr0_sq)_nw$(num_walkers).pdf"))
    end

    path = result_path(stage_dir, "num_walkers_convergence", (N=num_part, nr0sq=nr0_sq, num_walkers=num_walkers))
    save_run(path, (; num_part, nr0_sq, L, R_opt, Δτ_chosen, num_walkers, num_steps, num_equil, τ_total, quadratic),
             result; overwrite=force_rerun_overwrite)
end

x = 1.0 ./ num_walkers_vals
intercept, intercept_err, coeffs, fitted_fn = weighted_extrapolation(x, E_vals, E_err_vals; powers=[1,2])
println("\nExtrapolated E (1/N_walkers → 0) = $intercept ± $intercept_err")

x_range = collect(LinRange(0, maximum(x), 100))
p = scatter(x, E_vals; yerror=E_err_vals, xlabel=L"1/N_{\mathrm{walkers}}", ylabel=L"E",
            title="Walker-count convergence, N=$(num_part), nr0²=$(nr0_sq)", label="DMC data", framestyle=:box)
plot!(x_range, fitted_fn.(x_range); label="weighted fit", linestyle=:dash)
scatter!([0.0], [intercept]; yerror=[intercept_err], label="extrapolated (1/N→0)", marker=:star5, markersize=10, color=:red)
display(p)
mkpath(joinpath(stage_dir, "data", "plots"))
savefig(p, joinpath(stage_dir, "data", "plots", "num_walkers_convergence_N$(num_part)_nr0sq$(nr0_sq).pdf"))

final_path = result_path(stage_dir, "DMC", (N=num_part, nr0sq=nr0_sq))
save_run(final_path, (; num_part, nr0_sq, L, R_opt, Δτ_chosen, num_steps, num_equil, τ_total, quadratic),
         (; E_extrapolated=intercept, E_extrapolated_err=intercept_err, num_walkers_vals, E_vals, E_err_vals);
         overwrite=force_rerun_overwrite)
println("Saved final DMC result to: ", final_path)