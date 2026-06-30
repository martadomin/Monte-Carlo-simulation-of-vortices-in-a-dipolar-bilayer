"""
find_delta_bilayer.jl

Standalone script to find the optimal delta for move_all VMC in the bilayer system.
Reads R0_opt, R_match, energy_b from the R0 sweep file produced by optimize_R0.jl.
Does NOT call optimize_R0.jl or run_VMC_bilayer.jl.

Sets on exit:
  delta_opt_moveall, Δτ_DMC
"""

using Random, DelimitedFiles, Plots, LaTeXStrings, ProgressMeter

include(normpath(joinpath(@__DIR__, "src", "utils.jl")))
include(normpath(joinpath(@__DIR__, "src", "jastrow.jl")))
include(normpath(joinpath(@__DIR__, "src", "energy.jl")))
include(normpath(joinpath(@__DIR__, "src", "metropolis.jl")))
include(normpath(joinpath(@__DIR__, "src", "shooting_method.jl")))
include(normpath(joinpath(@__DIR__, "src", "observables.jl")))

# ── Parameters ─────────────────────────────────────────────────────────────────
N      = 60
nr0sq  = 1.0
h      = 0.3
L      = sqrt(N / nr0sq)

r_min     = 1e-6
Δ_shoot   = 1e-4
tol_shoot = 1e-10

num_sweep_steps = 10^4
delta_vals = [0.01, 0.02, 0.05, 0.08, 0.10, 0.12, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50]

# ── Read R0_opt, R_match, energy_b from sweep file ─────────────────────────────
sweep_path = joinpath(@__DIR__, "data", "sweep_results",
             "R0_sweep_Stage2_N$(N)_nr0sq$(nr0sq)_h$(h).txt")

# R0_opt   = NaN
# R_match  = NaN
# energy_b = NaN

# open(sweep_path, "r") do io
#     for line in eachline(io)
#         # strip whitespace so "R0_opt      = ..." and "R0_opt= ..." both work
#         s = strip(line)
#         if startswith(s, "R0_opt")
#             R0_opt   = parse(Float64, strip(split(s, "=")[2]))
#         elseif startswith(s, "energy_b")
#             energy_b = parse(Float64, strip(split(s, "=")[2]))
#         elseif startswith(s, "#") && occursin("R_match", s)
#             m = match(r"R_match=([\d.eE+\-]+)", s)
#             m !== nothing && (R_match = parse(Float64, m.captures[1]))
#         end
#     end
# end

R0_opt = 0.7281208690869945
R_match = 0.8349980438651612


@assert !isnan(R0_opt)  "R0_opt not found in $sweep_path"
@assert !isnan(R_match) "R_match not found in $sweep_path"

println("="^70)
println("DELTA SWEEP (move_all)  h = $h r₀")
println("="^70)
println("  R0_opt   = $R0_opt")
println("  R_match  = $R_match")
println("  L        = $L\n")

Constants = calculate_constants(L, R_match)

# ── Build fAB ─────────────────────────────────────────────────────────────────
_, _, itp_u, itp_up, itp_upp, _ =
    build_fAB(h, R0_opt, r_min, Δ_shoot, tol_shoot, nr0sq, N)

# ── Initial configuration ──────────────────────────────────────────────────────
x_coord, y_coord = random_initial_config(N, L, "Uniform")
x_A = x_coord[1:N÷2];     y_A = y_coord[1:N÷2]
x_B = x_coord[N÷2+1:end]; y_B = y_coord[N÷2+1:end]

# ── Sweep ──────────────────────────────────────────────────────────────────────
acceptance_vals = Float64[]

for δ in delta_vals
    _, _, _, _, _, _, _, _, _, acc, _, _ = metropolis(
        N, num_sweep_steps, δ, L, h, R_match, R0_opt,
        itp_u, itp_up, itp_upp, Constants;
        x_A_init = x_A, y_A_init = y_A,
        x_B_init = x_B, y_B_init = y_B,
        progress = true,
        move_all = true
    )
    push!(acceptance_vals, acc)
    println("  δ = $(rpad(string(δ), 6)) → acceptance = $(round(acc * 100, digits=2))%")
end

# ── Find optimal delta ─────────────────────────────────────────────────────────
best_idx          = argmin(abs.(acceptance_vals .- 0.5))
delta_opt_moveall = delta_vals[best_idx]
acc_opt           = acceptance_vals[best_idx]
Δτ_DMC            = delta_opt_moveall^2   # D = 1/2, Δτ = δ²/(2D) = δ²

println("\n  Optimal δ = $delta_opt_moveall  (acceptance = $(round(acc_opt * 100, digits=2))%)")
println("  Δτ_DMC    = $Δτ_DMC")

# ── Plot ───────────────────────────────────────────────────────────────────────
p = plot(delta_vals, acceptance_vals .* 100;
         xlabel  = L"\delta",
         ylabel  = "Acceptance (%)",
         title   = "move_all delta sweep  N=$N, h=$h r₀",
         marker  = :circle,
         xscale  = :log10,
         label   = "acceptance",
         legend  = :topright)
hline!(p, [50.0]; linestyle = :dash, color = :red,   label = "50%")
vline!(p, [delta_opt_moveall]; linestyle = :dot, color = :black, label = "optimal δ")
display(p)

# ── Save ───────────────────────────────────────────────────────────────────────
out_path = joinpath(@__DIR__, "data", "sweep_results",
           "delta_moveall_N$(N)_nr0sq$(nr0sq)_h$(h).txt")
open(out_path, "w") do io
    println(io, "# move_all delta sweep — bilayer")
    println(io, "# N=$(N), nr0sq=$(nr0sq), h=$(h), R0_opt=$(R0_opt), R_match=$(R_match)")
    println(io, "# delta_opt=$(delta_opt_moveall), Δτ_DMC=$(Δτ_DMC)")
    println(io, "#")
    println(io, "# delta\tacceptance")
    for (δ, acc) in zip(delta_vals, acceptance_vals)
        println(io, "$(δ)\t$(acc)")
    end
    println(io, "#")
    println(io, "# delta_opt\tDelta_tau_DMC")
    println(io, "$(delta_opt_moveall)\t$(Δτ_DMC)")
end
println("✓ Saved to: $out_path")