using DelimitedFiles

# ── Load exact binding energies ────────────────────────────────────
exact_data = readdlm(normpath(joinpath(@__DIR__, "..", "data", "binding_energy_dimer",
                     "dimer_binding_energy.txt")), '\t', skipstart=1)
h_exact  = Float64.(exact_data[:, 1])
eb_exact = Float64.(exact_data[:, 2])

# ── Load VMC sweep results ─────────────────────────────────────────
h_found        = Float64[]
E_per_N_vals   = Float64[]
E_minus_eb_vals = Float64[]
err_vals       = Float64[]

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
        occursin("=", s) || continue
        k, v = strip.(split(s, "=", limit=2))
        k == "R0_opt"   && (R0_opt  = parse(Float64, v))
        k == "E_opt"    && (E_opt   = parse(Float64, v))
        k == "err_opt"  && (err_opt = parse(Float64, v))
        k == "energy_b" && (eps_b   = parse(Float64, v))
    end
    isnan(L) && (L = sqrt(N / nr0sq))
    return R0_opt, E_opt, err_opt, eps_b, L
end

for h in h_vals
    local sweep_path, R0_opt, E_opt, err_opt, eps_b, L_file
    local n_layer, rc, tail, eb_exact_h

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

    # interpolate exact ε_b at this h
    i       = searchsortedlast(h_exact, h)
    i       = clamp(i, 1, length(h_exact) - 1)
    t       = (h - h_exact[i]) / (h_exact[i+1] - h_exact[i])
    eb_exact_h = (1-t) * eb_exact[i] + t * eb_exact[i+1]

    E_corr = E_opt

    push!(h_found,         h)
    push!(E_per_N_vals,    E_corr)
    push!(E_minus_eb_vals, E_corr - eb_exact_h / 2)
    push!(err_vals,        err_opt)
end

# ── Plot ───────────────────────────────────────────────────────────
pgfplotsx()

p = scatter(
    h_found, E_minus_eb_vals;
    yerror        = err_vals,
    xlabel        = L"h/r_0",
    ylabel        = L"E/N - \varepsilon_b/2 \, [\hbar^2/(mr_0^2)]",
    label         = L"VMC\ E/N",
    marker        = :circle,
    markersize    = 4,
    color         = :dodgerblue,
    grid          = true,
    gridalpha     = 0.25,
    framestyle    = :box,
    size          = (900, 600),
    tickfontsize  = 11,
    guidefontsize = 13,
    xlims         = (0.0, 1.0),
)

display(p)
savefig(p, normpath(joinpath(@__DIR__, "..", "data", "plots",
        "E_per_N_and_binding_energy_N$(N)_nr0sq$(nr0sq).pdf")))
readline()