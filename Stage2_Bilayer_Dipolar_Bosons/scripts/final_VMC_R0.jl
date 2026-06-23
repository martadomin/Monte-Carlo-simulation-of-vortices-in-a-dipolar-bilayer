# scripts/plot_E_vs_h.jl
# Stage 2 bilayer VMC: E/N as a function of layer separation h, with tail correction.
#
# Reads the "--- Optimal ---" block of each R0 sweep file:
#   data/sweep_results/R0_sweep_Stage2_N$(N)_nr0sq$(nr0sq)_h$(h).txt
#
# Tail correction (in-plane cutoff r_c = L/2, g(r) -> n_layer beyond cutoff):
#   E_tail/N = 2*pi*n_l/L  +  pi*n_l*(L/2)^2 / ((L/2)^2 + h^2)^(3/2)
#   intralayer (AA+BB)        interlayer (AB), from
#   V_AB(r) = (r^2 - 2h^2)/(r^2 + h^2)^(5/2) integrated over [L/2, inf)

ENV["PATH"] = "C:\\Users\\marta\\AppData\\Local\\Programs\\MiKTeX\\miktex\\bin\\x64;" * ENV["PATH"]
using DelimitedFiles, Plots, LaTeXStrings, PGFPlotsX

# ── CONFIGURATION ──────────────────────────────────────────────────────────────

subtract_binding = false      # true: also plot E/N + E_tail - |eps_b|/2 (Macia Fig. 1 style)

file_for_h(h) = joinpath(@__DIR__, "..", "data", "sweep_results",
                "R0_sweep_Stage2_N$(num_part)_nr0sq$(nr0_sq)_h$(h).txt")
# ───────────────────────────────────────────────────────────────────────────────

"""
    E_tail_per_N(L, h, n_layer)

Tail correction per particle for the bilayer with in-plane cutoff r_c = L/2.
"""
function E_tail_per_N(L::Float64, h::Float64, n_layer::Float64)
    rc = L / 2
    return 2π * n_layer / L + π * n_layer * rc^2 / (rc^2 + h^2)^(3/2)
end

"""
    read_optimal(path)

Parse the `--- Optimal ---` key = value block and the L from the header comment
of an R0 sweep file. Returns (R0_opt, E_opt, err_opt, eps_b, L).
"""
function read_optimal(path::String)
    R0_opt = E_opt = err_opt = eps_b = L = NaN
    for line in eachline(path)
        s = strip(line)
        isempty(s) && continue
        if startswith(s, "#")
            # header comment: "# L=7.745966692414834, L/2=3.872983346207417"
            if isnan(L) && occursin("L=", s)
                for tok in split(strip(s, ['#', ' ']), ",")
                    kv = split(strip(tok), "=")
                    length(kv) == 2 && strip(kv[1]) == "L" &&
                        (L = parse(Float64, strip(kv[2])))
                end
            end
            continue
        end
        occursin("=", s) || continue
        k, v = strip.(split(s, "=", limit=2))
        k == "R0_opt"   && (R0_opt  = parse(Float64, v))
        k == "E_opt"    && (E_opt   = parse(Float64, v))
        k == "err_opt"  && (err_opt = parse(Float64, v))
        k == "energy_b" && (eps_b   = parse(Float64, v))
    end
    isnan(L) && (L = sqrt(N / nr0_sq))   # fallback
    return R0_opt, E_opt, err_opt, eps_b, L
end

# # ── Load data ──────────────────────────────────────────────────────────────────
# h_found   = Float64[]
# E_per_N   = Float64[]   # E_opt is already E/N in the sweep files
# err_per_N = Float64[]
# eps_b_vals  = Float64[]
# E_tail_vals = Float64[]
# R0_opt_vals = Float64[]

# for h in h_vals
#     path = file_for_h(h)
#     if !isfile(path)
#         println("Skipping h = $h — file not found: $path")
#         continue
#     end
#     R0_opt, E_opt, err_opt, eps_b, L = read_optimal(path)
#     if isnan(E_opt)
#         println("Skipping h = $h — no Optimal block found in $path")
#         continue
#     end

#     n_layer = (num_part / 2) / L^2

#     push!(h_found,     h)
#     push!(E_per_N,     E_opt)
#     push!(err_per_N,   err_opt)
#     push!(eps_b_vals,  eps_b)
#     push!(E_tail_vals, E_tail_per_N(L, h, n_layer))
#     push!(R0_opt_vals, R0_opt)
# end

# isempty(h_found) && error("No result files found — check file_for_h() and h_vals.")

# E_corr = E_per_N .- eps_b_vals ./ 2   # E/N + E_tail - |eps_b|/2

# println("\n h/r0    R0_opt     E/N (raw)      E_tail/N      E/N (corr)     eps_b")
# for i in eachindex(h_found)
#     println(" $(rpad(h_found[i], 6))  $(rpad(round(R0_opt_vals[i], digits=4), 8))  ",
#             "$(rpad(round(E_per_N[i], digits=6), 12))  ",
#             "$(rpad(round(E_tail_vals[i], digits=6), 12))  ",
#             "$(rpad(round(E_corr[i], digits=6), 12))  ",
#             "$(round(eps_b_vals[i], digits=4))")
# end

# # ── Plot ───────────────────────────────────────────────────────────────────────
# # pgfplotsx()
# gr()


# col_raw  = "#0072B2"
# col_corr = "#D55E00"
# col_bind = "#009E73"

# p = plot(
#     xlabel = L"h/r_0",
#     ylabel = L"E/N\ \left[\hbar^2/(m r_0^2)\right]",
#     title  = L"\mathrm{Bilayer\ VMC},\ N = %$(num_part),\ nr_0^2 = %$(nr0_sq)",
#     framestyle = :box,
#     legend = :bottomright,
#     grid = true,
#     gridalpha = 0.3,
#     xlims = (0.0, maximum(h_found) + 0.1)
# )

# # plot!(p, h_found, eps_b_vals,
# #          yerror = err_per_N,
# #          label  = L"\mathrm{VMC,\ raw}",
# #          marker = :circle,
# #          markersize = 0,
# #          linewidth = 15,
# #          color  = col_raw)

# # scatter!(p, h_found, E_per_N,
# #          yerror = err_per_N,
# #          label  = L"\mathrm{VMC,\ raw}",
# #          marker = :circle,
# #          markersize = 5,
# #          color  = col_raw)

# scatter!(p, h_found, E_corr,
#          yerror = err_per_N,
#          label  = L"\mathrm{VMC} + E_\mathrm{tail}",
#          marker = :square,
#          markersize = 5,
#          color  = col_corr)

# # hline!(p, h_found, zeros(length(h_found)),
# #        label = L"\varepsilon_b\ \mathrm{(trial\ value\ at\ }R_0^\mathrm{opt}\mathrm{)}",
# #        linestyle = :dash,
# #        color = col_bind)
# # hline!(p, h_found, 4.6*ones(length(h_found)),
# #        label = L"\varepsilon_b\ \mathrm{(trial\ value\ at\ }R_0^\mathrm{opt}\mathrm{)}",
# #        linestyle = :dash,
# #        color = col_bind)

# # if subtract_binding
# #     # Macia et al. 2014, Fig. 1 style: E/N - |eps_b|/2
# #     # NOTE: eps_b here is the trial-function value at R0_opt, not the exact
# #     # dimer binding energy — see caveat in the accompanying notes.
# #     scatter!(p, h_found, E_corr .- abs.(eps_b_vals) ./ 2,
# #              yerror = err_per_N,
# #              label  = L"\mathrm{VMC} + E_\mathrm{tail} - |\varepsilon_b|/2",
# #              marker = :diamond,
# #              markersize = 5,
# #              color  = col_bind)
# # end

# plots_dir = joinpath(@__DIR__, "..", "data", "plots")
# isdir(plots_dir) || mkpath(plots_dir)
# savefig(p, joinpath(plots_dir, "E_vs_h_N$(num_part)_nr0sq$(nr0_sq).pdf"))
# savefig(p, joinpath(plots_dir, "E_vs_h_N$(num_part)_nr0sq$(nr0_sq).png"))
# println("\nSaved plot: E_vs_h_N$(num_part)_nr0sq$(nr0_sq).pdf")
# display(p)
# readline()