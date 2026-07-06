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
h_found        = Float64[]
E_raw_vals     = Float64[]   # no tail correction
E_corr_vals    = Float64[]   # with tail correction
err_vals       = Float64[]

for h in h_vals
    sweep_path = joinpath(@__DIR__, "..", "data", "sweep_results",
                          "R0_sweep_Stage2_N$(N)_nr0sq$(nr0sq)_h$(h).txt")
    if !isfile(sweep_path)
        println("Warning: file not found for h = $h")
        continue
    end
    R0_opt, E_opt, err_opt, eps_b, L_file = read_optimal(sweep_path)
    isnan(E_opt) && continue

    E_tail = tail_energy(nr0sq, N, h)

    println("h = $(rpad(h,5))  E/N = $(round(E_opt,digits=4))  " *
            "E_tail = $(round(E_tail,digits=4))  " *
            "E/N+tail = $(round(E_opt+E_tail,digits=4))")

    push!(h_found,     h)
    push!(E_raw_vals,  E_opt)
    push!(E_corr_vals, E_opt + E_tail)
    push!(err_vals,    err_opt)
end

# ── Single-layer reference ─────────────────────────────────────────
sweep_path_single_layer = joinpath(@__DIR__, "..", "..", "Stage1_2D_Dipol_System",
                                   "data", "sweep_results",
                                   "Rmatch_sweep_N$(N÷2)_nr0sq$(nr0sq/2).txt")

R_opt_single, E_opt_single_total = read_optimal_stage1(sweep_path_single_layer)

E_single_raw  = NaN
E_single_corr = NaN

if !isnan(E_opt_single_total)
    E_single_raw  = E_opt_single_total / (N ÷ 2)
    E_tail_single = 2π * (nr0sq/2)^(3/2) / sqrt(N ÷ 2)
    E_single_corr = E_single_raw + E_tail_single
    println("Single-layer E/N (raw)  = $(round(E_single_raw,  digits=4))")
    println("Single-layer E/N (corr) = $(round(E_single_corr, digits=4))")
else
    println("Warning: could not parse single-layer optimum from $sweep_path_single_layer")
end

# ── Plot: ────────────────────────────────────
pgfplotsx()

p = scatter(h_found, E_corr_vals;
    yerror        = err_vals,
    xlabel        = L"h/r_0",
    ylabel        = L"E/N \, [\hbar^2/(mr_0^2)]",
    label         = L"VMC\ $E/N + E_{\mathrm{tail}}$",
    marker        = :circle,
    markersize    = 8,
    color         = :dodgerblue,
    grid          = true,
    gridalpha     = 0.25,
    framestyle    = :box,
    tickfontsize  = 15,
    guidefontsize = 17,
    legendfontsize = 15,
    legendposition = :topleft,
    xlims = (0.0, 1.51),
    ylims = (-2.01, 5.0),
    size  = (900, 600),
)

plot!(p, h_exact, eb_exact ./ 2;
    label = L"\varepsilon_b/2", linewidth = 2, color = :black)

hline!(p, [0.0];
    color = :black, linestyle = :dash, linewidth = 1.0, label = "")

!isnan(E_single_corr) && hline!(p, [E_single_corr];
    color = :red, linestyle = :dash, linewidth = 1.0,
    label = L"E/N\ nr_0^2 = 0.5\ \mathrm{(single\ layer)}")

display(p)

plots_dir = normpath(joinpath(@__DIR__, "..", "data", "plots"))
mkpath(plots_dir)
savefig(p, joinpath(plots_dir,
        "Inset_Fig1_N$(N)_nr0sq$(nr0sq).pdf"))
println("✓ Saved to data/plots/Inset_Fig1_N$(N)_nr0sq$(nr0sq).pdf")