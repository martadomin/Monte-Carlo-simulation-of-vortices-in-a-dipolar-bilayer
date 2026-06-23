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
    isnan(L) && (L = sqrt(N / nr0sq))   # fallback
    return R0_opt, E_opt, err_opt, eps_b, L
end

"""
    read_optimal_stage1(path)

Parse the trailing "# Optimal R_match" block of a Stage I Rmatch sweep file.
That block has the form:
    # Optimal R_match
    R_opt	E_opt
    <value>	<value>
i.e. a header row of *column names* (no '='), followed by one data row.
The stored E_opt is the TOTAL energy (not per particle).
Returns (R_opt, E_opt_total) as Float64, or (NaN, NaN) if not found.
"""
function read_optimal_stage1(path::String)
    lines = readlines(path)
    idx = findfirst(l -> occursin("Optimal R_match", l), lines)
    idx === nothing && return NaN, NaN
    # idx -> "# Optimal R_match", idx+1 -> "R_opt\tE_opt", idx+2 -> values
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
h_found      = Float64[]
E_per_N_vals = Float64[]
err_vals     = Float64[]

for h in h_vals
    sweep_path = joinpath(@__DIR__, "..", "data", "sweep_results",
                          "R0_sweep_Stage2_N$(N)_nr0sq$(nr0sq)_h$(h).txt")
    if !isfile(sweep_path)
        println("Warning: file not found for h = $h")
        continue
    end
    R0_opt, E_opt, err_opt, eps_b, L_file = read_optimal(sweep_path)
    isnan(E_opt) && continue

    n_layer = (N/2) / L_file^2
    rc      = L_file / 2

    push!(h_found,      h)
    push!(E_per_N_vals, E_opt)
    push!(err_vals,     err_opt)
end

# ── Plot ───────────────────────────────────────────────────────────
pgfplotsx()

p = scatter(
    h_found, E_per_N_vals;
    yerror     = err_vals,
    xlabel     = L"h/r_0",
    ylabel     = L"E/N \, [\hbar^2/(mr_0^2)]",
    label      = L"VMC\ E/N",
    marker     = :circle,
    markersize = 4,
    color      = :dodgerblue,
    grid       = true,
    gridalpha  = 0.25,
    framestyle = :box,
    size       = (900, 600),
    tickfontsize  = 11,
    guidefontsize = 13,
    legendfontsize = 11,
    legendposition = :topleft,
    xlims = (0.0, 1.01),
    ylims = (-2.01, 4.5),
)

plot!(p,
    h_exact, eb_exact ./ 2;
    label     = L"\varepsilon_b/2",
    linewidth = 2,
    linestyle = :solid,
    color     = :black,
)

hline!(p, [0.0], color = :black, linestyle = :dash, linewidth = 1., label = "")

# Plot the energy per particle for a single layer of N/2 particles
# (independent-layer limit, h → ∞), read from the Stage I Rmatch sweep file.
sweep_path_single_layer = joinpath(@__DIR__, "..", "..", "Stage1_2D_Dipol_System",
                                   "data", "sweep_results",
                                   "Rmatch_sweep_N$(N÷2)_nr0sq$(nr0sq/2).txt")

R_opt_single, E_opt_single_total = read_optimal_stage1(sweep_path_single_layer)

if isnan(E_opt_single_total)
    println("Warning: could not parse single-layer optimum from $sweep_path_single_layer")
else
    E_opt_single = E_opt_single_total / (N ÷ 2)   # file stores TOTAL energy → convert to E/N
    println("Single-layer E/N = $E_opt_single  (R_match = $R_opt_single)")
    hline!(p, [E_opt_single], color = :red, linestyle = :dash, linewidth = 1.,
           label = L"E/N\ (single\ layer)")
end

display(p)
savefig(p, normpath(joinpath(@__DIR__, "..", "data", "plots",
        "Inset_Fig1_N$(N)_nr0sq$(nr0sq).pdf")))