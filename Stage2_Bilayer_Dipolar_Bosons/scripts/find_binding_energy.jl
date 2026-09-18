# Stage2_Bilayer_Dipolar_Bosons/scripts/find_binding_energy.jl
#
# Expects: h_vals, r_min, Δ_shoot, tol_shoot, L, stage_dir

h_found, eb_found = Float64[], Float64[]
for h in h_vals
    R_inf = L / 2
    Δ_h = min(Δ_shoot, h^(3/2) / 20.0)
    eb = find_energy_b(r_min, h, R_inf, Δ_h, tol_shoot)
    push!(h_found, h); push!(eb_found, eb)
    println("h = $(round(h, digits=4)) → ε_b = $(round(eb, digits=6))")
end

path = result_path(stage_dir, "dimer_binding_energy", (L=L,))
save_run(path, (; h_vals, r_min, Δ_shoot, tol_shoot, L), (; h_vals=h_found, eb_vals=eb_found))
println("Saved to: ", path)