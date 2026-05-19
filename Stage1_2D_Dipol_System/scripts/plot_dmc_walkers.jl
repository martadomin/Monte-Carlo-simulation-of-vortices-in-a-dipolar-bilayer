# Plot DMC energy vs number of walkers for a given nr0sq
using DelimitedFiles, Plots, Statistics

# ----------------- User parameters -----------------
# Edit these to match the run you want to inspect. `nr0sq_str` must match
# the exact string used in the results filenames (e.g. "32.0" or "0.25").
num_part = 30
nr0sq_str = "16.0"         # string as used in filenames
type_dmc = "quadratic"     # "quadratic" or "linear"
# If empty, script will autodetect available walker files for the chosen nr0sq
num_walkers_list = Int[]
# If provided, will pick the row with closest Δτ to this value; set to `nothing`
# to use the first row in each file.
Δτ_target = nothing

# ----------------- Paths -----------------
data_dir = joinpath(@__DIR__, "..", "data", "results", "DMC")
plots_dir = joinpath(@__DIR__, "..", "data", "plots", "DMC")
mkpath(plots_dir)

# ----------------- Discover files -----------------
files = readdir(data_dir)
pattern = "^dmc_N" * string(num_part) * "_nr0sq" * nr0sq_str * "_Nw(\\d+)_" * type_dmc * "\\.txt" * "\$"
regex = Regex(pattern)
found = Dict{Int,String}()
for f in files
    m = match(regex, f)
    if m !== nothing
        nwalk = parse(Int, m.captures[1])
        found[nwalk] = joinpath(data_dir, f)
    end
end

if length(found) == 0
    error("No result files found for num_part=$(num_part), nr0sq=$(nr0sq_str), type=$(type_dmc) in $data_dir")
end

if length(num_walkers_list) == 0
    num_walkers_list = sort(collect(keys(found)))
else
    # filter to available files
    num_walkers_list = [n for n in num_walkers_list if haskey(found, n)]
    if length(num_walkers_list) == 0
        error("None of the requested walker counts were found for nr0sq=$(nr0sq_str)")
    end
end

# ----------------- Read data -----------------
walkers = Int[]
energies = Float64[]    # E per particle
errors = Float64[]
Δτs = Float64[]
nr0_sq_32 = num_part * parse(Float64, nr0sq_str)^(3/2)  # for normalisation

for n in num_walkers_list
    path = found[n]
    tbl = nothing
    try
        tbl = readdlm(path, '\t', Float64, skipstart=1)
    catch err
        @warn "Failed to read $path: $err"
        continue
    end
    if tbl === nothing || size(tbl, 1) == 0
        @warn "$path contains no data rows"
        continue
    end

    # tbl columns: num_walkers, Δτ, E_dmc, Error_E_dmc
    rows = eachrow(tbl)
    chosen = nothing
    if Δτ_target === nothing
        chosen = first(rows)
    else
        # pick row with closest Δτ
        best = nothing; bestd = Inf
        for r in rows
            d = abs(r[2] - Δτ_target)
            if d < bestd
                bestd = d; best = r
            end
        end
        chosen = best
    end

    push!(walkers, n)
    E = chosen[3]
    σ = chosen[4]
    push!(energies, E / nr0_sq_32)  # normalise per particle
    push!(errors, σ / nr0_sq_32)
    push!(Δτs, chosen[2])
end

# ----------------- Plot -----------------
if length(walkers) == 0
    error("No valid data rows found to plot.")
end

# Sort by walker count
order = sortperm(walkers)
walkers = walkers[order]
energies = energies[order]
errors = errors[order]
Δτs = Δτs[order]

p = plot!(1 ./walkers, energies; yerror=errors, xlabel="1 / Number of walkers",
            ylabel="E / N", title="DMC energy vs walkers (nr0sq=$(nr0sq_str))",
            marker=:circle, legend=false, xscale=:log10)

# Annotate Δτ used (if consistent across files, display it)
unique_Δτ = unique(Δτs)
if length(unique_Δτ) == 1
    annot = "Δτ=$(round(unique_Δτ[1], sigdigits=4))"
else
    annot = "Δτs vary: $(join(round.(unique_Δτ, sigdigits=4), ", "))"
end
annotate!(maximum(walkers), minimum(energies), text(annot, :right, 8))

out_path = joinpath(plots_dir, "dmc_walkers_nr0sq$(nr0sq_str).png")
savefig(p, out_path)
println("Saved plot to: $out_path")

# Print table
println("\nSummary (num_walkers, Δτ, E/N, Error/N):")
for (n, dt, e, er) in zip(walkers, Δτs, energies, errors)
    println("$n\t$(round(dt, sigdigits=6))\t$(round(e, digits=6))\t$(round(er, digits=6))")
end
