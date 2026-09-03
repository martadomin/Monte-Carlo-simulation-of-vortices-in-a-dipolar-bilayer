# main_DMC_bilayer.jl
using DelimitedFiles, Plots, LaTeXStrings, Statistics, ProgressMeter, Base.Threads

include(normpath(joinpath(@__DIR__, "src", "utils.jl")))
include(normpath(joinpath(@__DIR__, "src", "jastrow.jl")))
include(normpath(joinpath(@__DIR__, "src", "energy.jl")))
include(normpath(joinpath(@__DIR__, "src", "metropolis.jl")))
include(normpath(joinpath(@__DIR__, "src", "dmc.jl")))
include(normpath(joinpath(@__DIR__, "src", "shooting_method.jl")))
include(normpath(joinpath(@__DIR__, "src", "observables.jl")))

# ── Parameters ────────────────────────────────────────────────────────────────
N             = 60
nr0sq         = 1.0
num_walkers   = 500
quadratic     = true
type_dmc      = quadratic ? "quadratic" : "linear"
total_time    = 20.0
equil_time    = 4.0
τ_factor      = 0.1

h_vals_to_run = if !isempty(ARGS)
    [parse(Float64, ARGS[1])]
else
    collect(range(0.3, 1.5, step = 0.1))
end

r_min     = 1e-6
Δ_shoot   = 1e-4
tol_shoot = 1e-10

# # ── R_match (Stage I, AA/BB — h-independent) ────────────────────────────────
# vmc_path_stage1 = joinpath(@__DIR__, "..", "Stage1_2D_Dipol_System", "data", "results", "VMC",
#            "vmc_N$(N ÷ 2)_nr0sq$(nr0sq / 2).txt")
# vmc_data_stage1 = readdlm(vmc_path_stage1, '\t', String, skipstart=1)
# R_match = parse(Float64, vmc_data_stage1[1, 4])

# L = sqrt(N / nr0sq)

# # ── Main loop over h ─────────────────────────────────────────────────────────
# for h in h_vals_to_run

#     println("\n" * "="^70)
#     println("BILAYER DMC PRODUCTION RUN   h = $h r₀")
#     println("="^70)

#     # ── R0_opt + E_ref_initial (Stage II VMC results, production subfolder) ──
#     vmc_path_stage2 = joinpath(@__DIR__, "data", "results", "VMC", "production",
#                "VMC_Stage2_N$(N)_nr0sq$(nr0sq)_h$(h).txt")
#     vmc_data_stage2 = readdlm(vmc_path_stage2, '\t', String, skipstart=1)
#     R0_opt        = parse(Float64, vmc_data_stage2[1, 4])
#     E_ref_initial = parse(Float64, vmc_data_stage2[1, 5]) * N

#     Constants = calculate_constants(L, R_match)

#     # ── Δτ_ref from the move_all delta sweep (find_delta_bilayer.jl output) ─
#     delta_path = joinpath(@__DIR__, "data", "sweep_results",
#                  "delta_moveall_N$(N)_nr0sq$(nr0sq)_h$(h).txt")
#     delta_data = readdlm(delta_path, '\t', Float64; comments=true, comment_char='#')
#     Δτ_ref = delta_data[end, 2]   # last row = [delta_opt  Δτ_DMC] summary line

#     # ── Build fAB (cubic-spline interlayer Jastrow derivatives) ─────────────
#     _, _, _, itp_up, itp_upp, _ =
#         build_fAB(h, R0_opt, r_min, Δ_shoot, tol_shoot, nr0sq, N)

#     # ── Load VMC equilibrated configuration (A/B split) ──────────────────────
#     config_path = joinpath(@__DIR__, "data", "results", "VMC", "configs",
#                   "VMC_final_config_N$(N)_nr0sq$(nr0sq)_h$(h).txt")

#     config_lines = readlines(config_path)
#     layer_col = String[]
#     x_coord   = Float64[]
#     y_coord   = Float64[]
#     for line in config_lines
#         s = strip(line)
#         isempty(s) && continue
#         startswith(s, "#") && continue
#         cols = split(s, '\t')
#         length(cols) == 3 || continue
#         push!(layer_col, cols[1])
#         push!(x_coord, parse(Float64, cols[2]))
#         push!(y_coord, parse(Float64, cols[3]))
#     end

#     @assert length(x_coord) == N "Expected $N particles in $config_path, found $(length(x_coord))"
#     @assert all(==( "A"), layer_col[1:N÷2])     "Layer A block is not contiguous at the start"
#     @assert all(==( "B"), layer_col[N÷2+1:end]) "Layer B block is not contiguous at the end"

#     N_half = N ÷ 2
#     xA_init = x_coord[1:N_half];     yA_init = y_coord[1:N_half]
#     xB_init = x_coord[N_half+1:end]; yB_init = y_coord[N_half+1:end]

#     println("  R0_opt        = $R0_opt")
#     println("  R_match       = $R_match")
#     println("  E_ref_initial = $(E_ref_initial/N)")
#     println("  Δτ_ref        = $Δτ_ref")

#     # ── DMC run (fixed τ, fixed walker count) ────────────────────────────────
#     results_path = joinpath(@__DIR__, "data", "results", "DMC",
#                    "DMC_Stage2_N$(N)_nr0sq$(nr0sq)_h$(h).txt")

#     Δτ        = Δτ_ref * τ_factor
#     num_steps = max(10^4, round(Int, total_time / Δτ))
#     num_equil = max(2000, round(Int, equil_time / Δτ))

#     println("\n═══ h=$h | Δτ=$(round(Δτ, sigdigits=3)) | steps=$num_steps | walkers=$num_walkers ═══")

#     E_dmc, E_dmc_err, E_history = dmc(
#         xA_init, yA_init, xB_init, yB_init,
#         num_walkers, N, num_steps,
#         Δτ, L, h, R_match, Constants, R0_opt, itp_up, itp_upp,
#         E_ref_initial, num_walkers;
#         num_equil   = num_equil,
#         quadratic   = quadratic,
#         plot_energy = false
#     )

#     # Block averaging
#     block_sizes = [10, 20, 30, 40, 50, 100, 150, 200,
#                    300, 400, 500, 600, 700, 800, 900, 1000,
#                    1100, 1200, 1300, 1400, 1500, 1600, 1700,
#                    1800, 1900, 2000]
#     sigmas = [blocking_statistics(E_history, B)[2] for B in block_sizes]

#     plateau = detect_plateau(block_sizes, sigmas, window_size=4, rtol=0.02)

#     avg_E, σ = blocking_statistics(E_history, plateau)

#     println("E/N = $(round(avg_E/N, digits=5)) ± $(round(σ/N, digits=5))")

#     open(results_path, "w") do io
#         println(io, "num_walkers\tΔτ\tE_dmc\tError_E_dmc")
#         println(io, "$(num_walkers)\t$(Δτ)\t$(avg_E)\t$(σ)")
#     end

#     println("Saved: $results_path")
# end

# Create relevant plots
include(normpath(joinpath(@__DIR__, "scripts", "plot_fig1.jl")))

include(normpath(joinpath(@__DIR__, "scripts", "plot_inset_fig1.jl")))

println("\nAll bilayer DMC production runs completed.")