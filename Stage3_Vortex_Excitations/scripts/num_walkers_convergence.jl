# Stage3_Vortex_Excitations/scripts/num_walkers_convergence.jl
#
# Expects from main_DMC.jl: num_part, nr0_sq, h, lA, lB, d, L, r_min,
# Δ_shoot, tol_shoot, Δτ_chosen, num_walkers_vals, num_steps_dmc,
# quadratic, stage_dir, stage1_dir, stage2_dir

using Plots, LaTeXStrings

N_half = num_part ÷ 2
R_match = load_Rmatch(stage1_dir, N_half, L)
Constants = calculate_constants(L, R_match)
R0 = load_R0(stage2_dir, num_part, nr0_sq, h)

fAB_result = build_fAB(h, R0, r_min, Δ_shoot, tol_shoot, nr0_sq, num_part)
fAB_result === nothing && error("build_fAB failed at R0=$R0 — cannot run walker-count convergence.")
_, _, itp_u, itp_up, itp_upp, _ = fAB_result

x_vortex_A, y_vortex_A = L/2, L/2
x_vortex_B, y_vortex_B = L/2 + d, L/2

trial = VortexBilayerTrial(; R_match, Constants, R0, itp_u, itp_up, itp_upp, h, lA, lB,
                              x_vortex_A, y_vortex_A, x_vortex_B, y_vortex_B)

vmc_run = load_run(result_path(stage_dir, "VMC_vortex", (N=num_part, nr0sq=nr0_sq, h=h, lA=lA, lB=lB, d=d)); run=1)
E_ref_initial = vmc_run.result.avg_energy

coords_init = init_random_config((A=N_half, B=N_half), L)

E_vals, E_err_vals = Float64[], Float64[]
best_observables = nothing   # captured from the largest num_walkers — closest to 1/N→0

for num_walkers in num_walkers_vals
    println("\n--- num_walkers = $num_walkers, d=$d ---")
    result = dmc(trial, coords_init, num_walkers, num_steps_dmc, Δτ_chosen, L,
                 E_ref_initial, num_walkers; num_equil=num_steps_dmc ÷ 5, quadratic=quadratic)
    push!(E_vals, result.E_dmc)
    push!(E_err_vals, result.E_dmc_err)
    println("E = $(result.E_dmc) ± $(result.E_dmc_err)")

    num_walkers == maximum(num_walkers_vals) && (best_observables = result.observables)

    if plot_dmc_trace_diagnostic
        p_trace = plot_dmc_trace(result, num_steps_dmc ÷ 5;
                                  title="num_walkers=$(num_walkers), N=$(num_part), nr0²=$(nr0_sq), h=$(h), d=$(round(d,digits=4))")
        display(p_trace)
        savefig(p_trace, joinpath(stage_dir, "data", "plots",
                "dmc_trace_N$(num_part)_nr0sq$(nr0_sq)_h$(h)_d$(round(d,digits=4))_nw$(num_walkers).pdf"))
    end

    path = result_path(stage_dir, "num_walkers_convergence", (N=num_part, nr0sq=nr0_sq, h=h, lA=lA, lB=lB, d=d, num_walkers=num_walkers))
    save_run(path, (; num_part, nr0_sq, h, lA, lB, d, L, R_match, R0, Δτ_chosen, num_walkers, num_steps_dmc, quadratic), result)
end

x = 1.0 ./ num_walkers_vals
intercept, intercept_err, coeffs, fitted_fn = weighted_extrapolation(x, E_vals, E_err_vals; powers=[1,2])
println("\nExtrapolated E (1/N_walkers → 0) = $intercept ± $intercept_err")

x_range = collect(LinRange(0, maximum(x), 100))
p = scatter(x, E_vals; yerror=E_err_vals, xlabel=L"1/N_{\mathrm{walkers}}", ylabel=L"E",
            title="Walker-count convergence, N=$(num_part), nr0²=$(nr0_sq), h=$(h), d=$(round(d,digits=4))",
            label="DMC data", framestyle=:box)
plot!(x_range, fitted_fn.(x_range); label="weighted fit", linestyle=:dash)
scatter!([0.0], [intercept]; yerror=[intercept_err], label="extrapolated (1/N→0)", marker=:star5, markersize=10, color=:red)
display(p)
savefig(p, joinpath(stage_dir, "data", "plots", "num_walkers_convergence_N$(num_part)_nr0sq$(nr0_sq)_h$(h)_d$(round(d,digits=4)).pdf"))

final_path = result_path(stage_dir, "DMC_vortex", (N=num_part, nr0sq=nr0_sq, h=h, lA=lA, lB=lB, d=d))
save_run(final_path, (; num_part, nr0_sq, h, lA, lB, d, L, R_match, R0, Δτ_chosen, num_steps_dmc, quadratic),
         (; E_extrapolated=intercept, E_extrapolated_err=intercept_err, num_walkers_vals, E_vals, E_err_vals,
            observables=best_observables))
println("Saved final DMC result to: ", final_path)