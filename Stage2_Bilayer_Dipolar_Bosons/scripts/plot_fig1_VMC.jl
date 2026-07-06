using DelimitedFiles

"""
    read_optimal(path)
Parse the `--- Optimal ---` key = value block and the L from the header comment
of a Stage II R0 sweep file. Returns (R0_opt, E_opt, err_opt, eps_b, L).
"""
function read_optimal(path::String)
    R0_opt = E_opt = err_opt = eps_b = L = NaN
    for line in eachline(path)
        s = strip(line)
        isempty(s) && continue
        if startswith(s, "#")
            if isnan(L) && occursin("L=", s)
                for tok in split(strip(s, ['#', ' ']), ",")
                    kv = split(strip(tok), "=")
                    length(kv) == 2 && strip(kv[1]) == "L" &&
                        (L = parse(Float64, strip(kv[2])))
                end
            end
            continue
        end
        occursin("=", s) && !occursin("\t", s) || continue
        k, v = strip.(split(s, "=", limit=2))
        k == "R0_opt"   && (R0_opt  = parse(Float64, v))
        k == "E_opt"    && (E_opt   = parse(Float64, v))
        k == "err_opt"  && (err_opt = parse(Float64, v))
        k == "energy_b" && (eps_b   = parse(Float64, v))
    end
    isnan(L) && (L = sqrt(N / nr0sq))
    return R0_opt, E_opt, err_opt, eps_b, L
end

"""
    read_optimal_stage1(path)
Parse the trailing "# Optimal R_match" block of a Stage I Rmatch sweep file.
Returns (R_opt, E_opt_total) as Float64, or (NaN, NaN) if not found.
"""
function read_optimal_stage1(path::String)
    lines = readlines(path)
    idx = findfirst(l -> occursin("Optimal R_match", l), lines)
    idx === nothing && return NaN, NaN
    length(lines) < idx + 2 && return NaN, NaN
    vals = split(strip(lines[idx + 2]), '\t')
    length(vals) < 2 && return NaN, NaN
    return parse(Float64, vals[1]), parse(Float64, vals[2])
end

# ── Load exact binding energies ────────────────────────────────────
exact_data = readdlm(normpath(joinpath(@__DIR__, "..", "data", "binding_energy_dimer",
                     "dimer_binding_energy.txt")), '\t', skipstart=1)
h_exact  = Float64.(exact_data[:, 1])
eb_exact = Float64.(exact_data[:, 2])

# ── Load VMC sweep results ─────────────────────────────────────────
h_found          = Float64[]
E_minus_eb_vals  = Float64[]   # (E/N + E_tail) - ε_b/2
err_vals         = Float64[]

for h in h_vals
    local sweep_path, R0_opt, E_opt, err_opt, eps_b, L_file
    local i, t, eb_exact_h, E_tail, E_corr

    sweep_path = joinpath(@__DIR__, "..", "data", "sweep_results",
                          "R0_sweep_Stage2_N$(N)_nr0sq$(nr0sq)_h$(h).txt")
    if !isfile(sweep_path)
        println("Warning: file not found for h = $h")
        continue
    end

    R0_opt, E_opt, err_opt, eps_b, L_file = read_optimal(sweep_path)
    isnan(E_opt) && continue

    # Tail correction
    E_tail = tail_energy(nr0sq, N, h)
    E_corr = E_opt + E_tail

    # Interpolate exact ε_b at this h
    i = searchsortedlast(h_exact, h)
    i = clamp(i, 1, length(h_exact) - 1)
    t = (h - h_exact[i]) / (h_exact[i+1] - h_exact[i])
    eb_exact_h = (1-t) * eb_exact[i] + t * eb_exact[i+1]

    println("h = $(rpad(h,5))  E/N = $(round(E_opt,digits=4))  " *
            "E_tail = $(round(E_tail,digits=4))  " *
            "E_corr = $(round(E_corr,digits=4))  " *
            "ε_b/2 = $(round(eb_exact_h/2,digits=4))  " *
            "E_corr - ε_b/2 = $(round(E_corr - eb_exact_h/2,digits=4))")

    push!(h_found,         h)
    push!(E_minus_eb_vals, E_corr - eb_exact_h / 2)
    push!(err_vals,        err_opt)
end

# ── Single-layer reference lines ──────────────────────────────────
# Macià uses two horizontal lines:
#   1. Single layer at nr0sq = 0.5 (h → ∞ limit, each layer decouples)
#   2. Single layer at nr0sq = 8*0.5 = 4.0 (tightly bound dimer limit, h → 0)
# Both are plotted as E/N_single - ε_b/2 where ε_b → 0 for large h

# Line 1: h → ∞ limit — single layer at nr0sq = nr0sq/2 = 0.5
sweep_path_inf = joinpath(@__DIR__, "..", "..", "Stage1_2D_Dipol_System",
                          "data", "sweep_results",
                          "Rmatch_sweep_N$(N÷2)_nr0sq$(nr0sq/2).txt")
R_opt_inf, E_opt_inf_total = read_optimal_stage1(sweep_path_inf)

E_ref_inf = NaN
if !isnan(E_opt_inf_total)
    E_tail_inf = 2π * (nr0sq/2)^(3/2) / sqrt(N÷2)
    E_ref_inf  = E_opt_inf_total / (N÷2) + E_tail_inf
    # subtract ε_b/2 → 0 for large h, so reference line is just E/N
    println("h→∞ reference (E/N + tail) = $(round(E_ref_inf, digits=4))")
else
    println("Warning: could not read h→∞ single-layer reference")
end

# Line 2: h → 0 limit — dimers behave as single layer of dipoles
# with effective dipole length r̃₀ = 8r₀ → nr̃₀² = nr₀² × 64
# but we only plot this if we have the data
sweep_path_dim = joinpath(@__DIR__, "..", "..", "Stage1_2D_Dipol_System",
                          "data", "sweep_results",
                          "Rmatch_sweep_N$(N÷2)_nr0sq$(nr0sq/2).txt")
R_opt_dim, E_opt_dim_total = read_optimal_stage1(sweep_path_dim)

E_ref_dim = NaN
if !isnan(E_opt_dim_total)
    E_tail_dim = 2π * (nr0sq/2)^(3/2) / sqrt(N÷2)
    E_ref_dim  = E_opt_dim_total / (N÷2) + E_tail_dim
    println("h→0  reference (dimer, E/N + tail) = $(round(E_ref_dim, digits=4))")
else
    println("Note: dimer limit reference file not found — skipping that line")
end

# ── Plot ───────────────────────────────────────────────────────────
pgfplotsx()

p = scatter(
    h_found, E_minus_eb_vals;
    yerror        = err_vals,
    xlabel        = L"h/r_0",
    ylabel        = L"E/N - \varepsilon_b/2 \, [\hbar^2/(mr_0^2)]",
    label = L"\mathrm{VMC}",
    marker        = :circle,
    markersize    = 8,
    color         = :dodgerblue,
    grid          = true,
    gridalpha     = 0.25,
    framestyle    = :box,
    size          = (900, 600),
    tickfontsize  = 15,
    guidefontsize = 17,
    legendfontsize = 15,
    legendposition = :topright,
    xlims = (0.0, 1.5),
    ylims = (3.5, 4.6)
)

# h → ∞ reference line
if !isnan(E_ref_inf)
    hline!(p, [E_ref_inf],
           color = :black, linestyle = :dash, linewidth = 1.5,
           label = L"E/N\ nr_0^2 = 0.5\ \mathrm{(single\ layer)}")
end

# # h → 0 dimer reference line
# if !isnan(E_ref_dim)
#     hline!(p, [E_ref_dim],
#            color = :gray, linestyle = :dash, linewidth = 1.5,
#            label = L"h \to 0\ \mathrm{limit\ (dimer)}")
# end

# # Zero line
# hline!(p, [0.0],
#        color = :black, linestyle = :dot, linewidth = 1., label = "")

display(p)

plots_dir = normpath(joinpath(@__DIR__, "..", "data", "plots"))
mkpath(plots_dir)
savefig(p, joinpath(plots_dir, "Fig1_N$(N)_nr0sq$(nr0sq).pdf"))
println("✓ Saved to data/plots/Fig1_N$(N)_nr0sq$(nr0sq).pdf")
readline()