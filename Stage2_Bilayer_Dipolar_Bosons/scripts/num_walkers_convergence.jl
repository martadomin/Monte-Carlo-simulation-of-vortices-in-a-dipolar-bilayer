# Stage2_Bilayer_Dipolar_Bosons/scripts/num_walkers_convergence.jl
#
# Expects from main_DMC.jl: num_part, nr0_sq, h, L, r_min, Δ_shoot, tol_shoot,
# Δτ_chosen (from dtau_convergence.jl), num_walkers_vals, num_steps_dmc,
# quadratic, stage_dir, stage1_dir

using Plots, LaTeXStrings

N_half = num_part ÷ 2
R_match = load_Rmatch(stage1_dir, N_half, L)
Constants = calculate_constants(L, R_match)

R0_run = load_run(result_path(stage_dir, "R0_optimum", (N=num_part, nr0sq=nr0_sq, h=h)); run=1)
R0_opt = R0_run.result.R0_opt

fAB_result = build_fAB(h, R0_opt, r_min, Δ_shoot, tol_shoot, nr0_sq, num_part)
fAB_result === nothing && error("build_fAB failed at R0_opt = $R0_opt — cannot run walker-count convergence.")
_, _, itp_u, itp_up, itp_upp, _ = fAB_result

trial = BilayerTrial(; R_match, Constants, R0=R0_opt, itp_u, itp_up, itp_upp, h)

vmc_run = load_run(result_path(stage_dir, "VMC", (N=num_part, nr0sq=nr0_sq, h=h)); run=1)
E_ref_initial = vmc_run.result.avg_energy

coords_init = init_random_config((A=N_half, B=N_half), L)

E_vals, E_err_vals = Float64[], Float64[]
for num_walkers in num_walkers_vals
    println("\n--- num_walkers = $num_walkers ---")
    result = dmc(trial, coords_init, num_walkers, num_steps_dmc, Δτ_chosen, L,
             E_ref_initial, num_walkers; num_equil=num_steps_dmc ÷ 5, quadratic=quadratic)
    push!(E_vals, result.E_dmc)
    push!(E_err_vals, result.E_dmc_err)
    println("E = $(result.E_dmc) ± $(result.E_dmc_err)")

    if plot_dmc_trace_diagnostic
        p_trace = plot_dmc_trace(result, num_steps_dmc ÷ 5; title="num_walkers=$(num_walkers), N=$(num_part), nr0²=$(nr0_sq), h=$(h)")
        display(p_trace)
        savefig(p_trace, joinpath(stage_dir, "data", "plots", "dmc_trace_N$(num_part)_nr0sq$(nr0_sq)_h$(h)_nw$(num_walkers).pdf"))
    end

    path = result_path(stage_dir, "num_walkers_convergence", (N=num_part, nr0sq=nr0_sq, h=h, num_walkers=num_walkers))
    save_run(path, (; num_part, nr0_sq, h, L, R_match, R0=R0_opt, Δτ_chosen, num_walkers, num_steps_dmc, quadratic), result)
end

x = 1.0 ./ num_walkers_vals
intercept, intercept_err, coeffs, fitted_fn = weighted_extrapolation(x, E_vals, E_err_vals; powers=[1,2])
println("\nExtrapolated E (1/N_walkers → 0) = $intercept ± $intercept_err")

x_range = collect(LinRange(0, maximum(x), 100))
p = scatter(x, E_vals; yerror=E_err_vals, xlabel=L"1/N_{\mathrm{walkers}}", ylabel=L"E",
            title="Walker-count convergence, N=$(num_part), nr0²=$(nr0_sq), h=$(h)", label="DMC data", framestyle=:box)
plot!(x_range, fitted_fn.(x_range); label="weighted fit", linestyle=:dash)
scatter!([0.0], [intercept]; yerror=[intercept_err], label="extrapolated (1/N→0)", marker=:star5, markersize=10, color=:red)
display(p)
savefig(p, joinpath(stage_dir, "data", "plots", "num_walkers_convergence_N$(num_part)_nr0sq$(nr0_sq)_h$(h).pdf"))

# Final result for this (nr0_sq, h) — this IS the production DMC answer.
final_path = result_path(stage_dir, "DMC", (N=num_part, nr0sq=nr0_sq, h=h))
save_run(final_path, (; num_part, nr0_sq, h, L, R_match, R0=R0_opt, Δτ_chosen, num_steps_dmc, quadratic),
         (; E_extrapolated=intercept, E_extrapolated_err=intercept_err, num_walkers_vals, E_vals, E_err_vals))
println("Saved final DMC result to: ", final_path)