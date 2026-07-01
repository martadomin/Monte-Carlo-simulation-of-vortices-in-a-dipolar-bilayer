# run_dmc_bilayer.jl
using DelimitedFiles, Plots, LaTeXStrings, Statistics, ProgressMeter, Base.Threads

include(normpath(joinpath(@__DIR__, "..", "src", "utils.jl")))
include(normpath(joinpath(@__DIR__, "..", "src", "jastrow.jl")))
include(normpath(joinpath(@__DIR__, "..", "src", "energy.jl")))
include(normpath(joinpath(@__DIR__, "..", "src", "metropolis.jl")))
include(normpath(joinpath(@__DIR__, "..", "src", "dmc.jl")))
include(normpath(joinpath(@__DIR__, "..", "src", "shooting_method.jl")))
include(normpath(joinpath(@__DIR__, "..", "src", "observables.jl")))

# ── Parameters ────────────────────────────────────────────────────────────────
N             = 60
nr0sq         = 1.0
h_vals        = range(0.3, 0.3, step = 0.1)           # sweep over interlayer separations
num_walkers_list = [100]
quadratic     = true
type_dmc      = quadratic ? "quadratic" : "linear"
total_time    = 20.0
equil_time    = 4.0
τ_factors     = [0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0]

r_min     = 1e-6
Δ_shoot   = 1e-4
tol_shoot = 1e-10

# ── R_match (Stage I, AA/BB — h-independent) ────────────────────────────────
# Path depth and filename convention must match find_delta_bilayer.jl exactly:
# two ".." to reach the sibling Stage1 project, N÷2 particles, "nr0sq" (no underscore).
vmc_path_stage1 = joinpath(@__DIR__, "..", "..", "Stage1_2D_Dipol_System", "data", "results", "VMC",
           "vmc_N$(N ÷ 2)_nr0sq$(nr0sq / 2).txt")
vmc_data_stage1 = readdlm(vmc_path_stage1, '\t', String, skipstart=1)
R_match = parse(Float64, vmc_data_stage1[1, 4])

L = sqrt(N / nr0sq)

# ── Main loop over h ─────────────────────────────────────────────────────────
for h in h_vals

    E_res = Float64[]
    E_err = Float64[]
    Δτ_res = Float64[]

    println("\n" * "="^70)
    println("BILAYER DMC CONVERGENCE STUDY   h = $h r₀")
    println("="^70)

    # ── R0_opt + E_ref_initial (Stage II VMC results, production subfolder) ──
    vmc_path_stage2 = joinpath(@__DIR__, "..", "data", "results", "VMC", "production",
               "VMC_Stage2_N$(N)_nr0sq$(nr0sq)_h$(h).txt")
    vmc_data_stage2 = readdlm(vmc_path_stage2, '\t', String, skipstart=1)
    R0_opt        = parse(Float64, vmc_data_stage2[1, 4])
    E_ref_initial = parse(Float64, vmc_data_stage2[1, 5]) * N

    Constants = calculate_constants(L, R_match)

    # ── Δτ_ref from the move_all delta sweep (find_delta_bilayer.jl output) ─
    delta_path = joinpath(@__DIR__, "..", "data", "sweep_results",
                 "delta_moveall_N$(N)_nr0sq$(nr0sq)_h$(h).txt")
    delta_data = readdlm(delta_path, '\t', Float64; comments=true, comment_char='#')
    Δτ_ref = delta_data[end, 2]   # last row = [delta_opt  Δτ_DMC] summary line

    # ── Build fAB (cubic-spline interlayer Jastrow derivatives) ─────────────
    _, _, _, itp_up, itp_upp, _ =
        build_fAB(h, R0_opt, r_min, Δ_shoot, tol_shoot, nr0sq, N)

    # ── Load VMC equilibrated configuration (A/B split) ──────────────────────
    # File format: "# ..." comment/metadata lines, then rows of
    # "layer\tx\ty" where layer ∈ {A, B}. Particles are already ordered
    # A-block then B-block, matching N_half split below.
    config_path = joinpath(@__DIR__, "..", "data", "results", "VMC", "configs",
                  "VMC_final_config_N$(N)_nr0sq$(nr0sq)_h$(h).txt")

    config_lines = readlines(config_path)
    layer_col = String[]
    x_coord   = Float64[]
    y_coord   = Float64[]
    for line in config_lines
        s = strip(line)
        isempty(s) && continue
        startswith(s, "#") && continue
        cols = split(s, '\t')
        length(cols) == 3 || continue
        push!(layer_col, cols[1])
        push!(x_coord, parse(Float64, cols[2]))
        push!(y_coord, parse(Float64, cols[3]))
    end

    @assert length(x_coord) == N "Expected $N particles in $config_path, found $(length(x_coord))"
    @assert all(==( "A"), layer_col[1:N÷2])     "Layer A block is not contiguous at the start"
    @assert all(==( "B"), layer_col[N÷2+1:end]) "Layer B block is not contiguous at the end"

    N_half = N ÷ 2
    xA_init = x_coord[1:N_half];     yA_init = y_coord[1:N_half]
    xB_init = x_coord[N_half+1:end]; yB_init = y_coord[N_half+1:end]

    println("  R0_opt        = $R0_opt")
    println("  R_match       = $R_match")
    println("  E_ref_initial = $(E_ref_initial/N)")
    println("  Δτ_ref        = $Δτ_ref")

    # ── Loop over requested walker counts ────────────────────────────────────
    for num_walkers in num_walkers_list
        results_path = joinpath(@__DIR__, "..", "data", "results", "DMC",
                       "convergence_dmc_bilayer_N$(N)_nr0sq$(nr0sq)_h$(h)_Nw$(num_walkers)_$(type_dmc).txt")

        open(results_path, "w") do io
            println(io, "num_walkers\tΔτ\tE_dmc\tError_E_dmc")

            for factor in τ_factors
                Δτ        = Δτ_ref * factor
                num_steps = max(10^4, round(Int, total_time / Δτ))
                num_equil = max(2000, round(Int, equil_time / Δτ))

                println("\n═══ h=$h | Δτ=$(round(Δτ, sigdigits=3)) | steps=$num_steps | walkers=$num_walkers ═══")

                E_dmc, E_dmc_err, E_history = dmc(
                    xA_init, yA_init, xB_init, yB_init,
                    num_walkers, N, num_steps,
                    Δτ, L, h, R_match, Constants, R0_opt, itp_up, itp_upp,
                    E_ref_initial, num_walkers;
                    num_equil   = num_equil,
                    quadratic   = quadratic,
                    plot_energy = false
                )

                # Block averaging
                block_sizes = [10, 20, 30, 40, 50, 100, 150, 200,
                               300, 400, 500, 600, 700, 800, 900, 1000,
                               1100, 1200, 1300, 1400, 1500, 1600, 1700,
                               1800, 1900, 2000]
                sigmas = [blocking_statistics(E_history, B)[2] for B in block_sizes]

                plateau = detect_plateau(block_sizes, sigmas, window_size=4, rtol=0.02)

                p = plot(block_sizes, sigmas, marker=:o,
                         xlabel="Block Size", ylabel="Std. Error",
                         title="Blocking Analysis (h=$h, Δτ=$(round(Δτ, sigdigits=3)))")
                vline!([plateau], label="Plateau Block Size = $plateau", linestyle=:dash, color=:red)
                display(p)

                avg_E, σ = blocking_statistics(E_history, plateau)

                println("E/N = $(round(avg_E/N, digits=5)) ± $(round(σ/N, digits=5))")

                push!(E_res, avg_E/N)
                push!(E_err, σ/N)
                push!(Δτ_res, Δτ)

                println(io, "$(num_walkers)\t$(Δτ)\t$(avg_E)\t$(σ)")
                flush(io)
            end

            # ── Plot convergence of E/N vs Δτ for this number of walkers ─────────
            sort_idx = sortperm(Δτ_res)
            p_conv = plot(Δτ_res[sort_idx], E_res[sort_idx]; yerror=E_err[sort_idx],
                xlabel = L"\Delta\tau",
                ylabel = L"E/N",
                title  = "DMC Convergence Study (h=$h, Nw=$num_walkers)",
                marker = :circle,
                linewidth = 2,
                label  = "E/N"
            )

            hline!(p_conv, [E_ref_initial / N]; linestyle=:dash, color=:red, label="VMC result")

            display(p_conv)

            readline()  # pause to view plot

            plot_path = joinpath(@__DIR__, "..", "data", "results", "DMC",
                        "convergence_plot_bilayer_N$(N)_nr0sq$(nr0sq)_h$(h)_Nw$(num_walkers)_$(type_dmc).pdf")

            savefig(p_conv, plot_path)

            println("Saved: $plot_path")

        end
        println("\nSaved: $results_path")
    end
end

println("\nAll bilayer DMC runs completed.")