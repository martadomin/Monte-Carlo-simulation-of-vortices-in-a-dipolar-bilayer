using DelimitedFiles, Plots, LaTeXStrings, PGFPlotsX, Interpolations

function RK4(f::Function, r0::Float64, y0::Vector{Float64}, Δ::Float64, h::Float64, energy_b::Float64)::Vector{Float64}
    k1 = f(r0, y0, h, energy_b)
    k2 = f(r0 + Δ/2, y0 + Δ/2 * k1, h, energy_b)
    k3 = f(r0 + Δ/2, y0 + Δ/2 * k2, h, energy_b)
    k4 = f(r0 + Δ, y0 + Δ * k3, h, energy_b)
    return y0 + (Δ/6) * (k1 + 2*k2 + 2*k3 + k4)
end

function V_AB(r::Float64, h::Float64)::Float64
    r2 = r^2
    h2 = h^2
    return (r2 - 2h2) / (r2 + h2)^(5/2)
end

function ode_rhs(r::Float64, y::Vector{Float64}, h::Float64, energy_b::Float64)::Vector{Float64}
    dy = zeros(2)
    dy[1] = y[2]
    dy[2] = -y[2]/r + (V_AB(r, h) - energy_b) * y[1]
    return dy
end

function initial_conditions(r0::Float64, h::Float64, energy_b::Float64)::Vector{Float64}
    VAB = V_AB(0.0, h)
    y1 = 1.0 + (VAB - energy_b) * r0^2 / 4
    y2 = (VAB - energy_b) * r0 / 2
    return [y1, y2]
end

function shooting_method(r0::Float64, h::Float64, energy_b::Float64, r_max::Float64, Δ::Float64)::Float64
    y0 = initial_conditions(r0, h, energy_b)
    r = r0
    while r < r_max
        y0 = RK4(ode_rhs, r, y0, Δ, h, energy_b)
        r += Δ
    end
    return y0[2]
end

function find_energy_b(r0::Float64, h::Float64, r_max::Float64,
                       Δ::Float64, tol::Float64)::Float64

    F(eb) = shooting_method(r0, h, eb, r_max, Δ)

    E_lo = -2.0 / h^3
    E_hi = -1e-14

    F_lo = F(E_lo)
    F_hi = F(E_hi)

    if F_lo * F_hi > 0
        error("No sign change in bracket (-2/h³, 0) for h=$h, r_max=$r_max. F_lo=$F_lo, F_hi=$F_hi")
    end

    while (E_hi - E_lo) > tol
        E_mid = (E_lo + E_hi) / 2.0
        F_mid = F(E_mid)
        if F_lo * F_mid < 0
            E_hi = E_mid
        else
            E_lo = E_mid
            F_lo = F_mid
        end
    end

    return (E_lo + E_hi) / 2.0
end

"""
    build_fAB(h, R0, r_min, Δ, tol, nr0sq, N)

Builds the interlayer Jastrow factor f_AB(r) for r ∈ [r_min, R0].

IMPORTANT: the dimer binding energy ε_b is a property of the isolated pair —
it must NOT depend on the variational cutoff R0. It is therefore always
computed with r_max = L/2 (half the simulation box), regardless of the R0
being swept. R0 only controls where the variational wavefunction is cut off
and renormalized; it does not change the underlying ε_b.
"""
function build_fAB(h::Float64, R0::Float64, r_min::Float64,
                   Δ::Float64, tol::Float64, nr0sq::Float64, N::Int)

    Δ = min(Δ, h^(3/2) / 20.0)
    L_box = sqrt(N / nr0sq)

    # ε_b variacional: eigenvalor que impone ψ'(R0) = 0
    # Distinto del ε_b físico excepto cuando R0 = L/2
    energy_b = find_energy_b(r_min, h, R0, Δ, tol)

    # ε_b físico: para reportar en output, usar r_max = L/2
    energy_b_phys = find_energy_b(r_min, h, L_box / 2, Δ, tol)

    n_steps = round(Int, (R0 - r_min) / Δ)
    r_grid  = Vector{Float64}(undef, n_steps + 1)
    psi     = Vector{Float64}(undef, n_steps + 1)
    dpsi    = Vector{Float64}(undef, n_steps + 1)

    y0        = initial_conditions(r_min, h, energy_b)
    r_grid[1] = r_min
    psi[1]    = y0[1]
    dpsi[1]   = y0[2]

    r = r_min
    for i in 1:n_steps
        y0          = RK4(ode_rhs, r, y0, Δ, h, energy_b)
        r          += Δ
        r_grid[i+1] = r
        psi[i+1]    = y0[1]
        dpsi[i+1]   = y0[2]
    end

    # Return nothing instead of erroring — caller decides what to do
    if psi[1] < 0.0 || any(i -> psi[i] * psi[i-1] < 0.0, 2:length(psi))
        return nothing
    end

    psi_R0 = psi[end]

    if abs(psi_R0) < 1e-6
        @warn "ψ(R0) is very close to zero (|ψ(R0)| = $(round(abs(psi_R0), sigdigits=3))) — " *
                "normalization will be unstable. Consider increasing R0 or adjusting Δ and tol."
        return nothing
    end

    psi  ./= psi_R0
    dpsi ./= psi_R0

    u_prime       = Vector{Float64}(undef, length(r_grid))
    u_doubleprime = Vector{Float64}(undef, length(r_grid))
    for i in 1:length(r_grid)
        u_prime[i]       = dpsi[i] / psi[i]
        d2psi_i          = (V_AB(r_grid[i], h) - energy_b) * psi[i] - dpsi[i]/r_grid[i]
        u_doubleprime[i] = d2psi_i / psi[i] - u_prime[i]^2
    end

    # Then build interpolants from these ODE-exact grid values
    r_range   = range(r_min, R0, length = length(r_grid))
    itp_u     = extrapolate(scale(interpolate(log.(psi), BSpline(Cubic(Natural(OnGrid())))), r_range), Flat())
    itp_up = extrapolate(scale(interpolate(u_prime, BSpline(Cubic(Flat(OnGrid())))), r_range), Flat())
    itp_upp   = extrapolate(scale(interpolate(u_doubleprime, BSpline(Cubic(Natural(OnGrid())))), r_range), Flat())

    return r_grid, psi, itp_u, itp_up, itp_upp, energy_b_phys
end

# # ════════════════════════════════════════════════════════════════════
# # SANITY CHECK — only runs when this file is executed directly
# # (julia .\src\shooting_method.jl), not when `include`d from elsewhere.
# # ════════════════════════════════════════════════════════════════════
# if abspath(PROGRAM_FILE) == @__FILE__

#     # ── Parameters (edit freely for ad-hoc checks) ─────────────────
#     r_min   = 1e-6
#     tol     = 1e-10
#     nr0sq   = 1.0
#     N       = 60
#     L_box   = sqrt(N / nr0sq)
#     R0_test = L_box / 2          # cutoff used for the f_AB build below

#     h_vals_check = range(0.3, 1.5, step = 0.1)   # range(start, stop; step=...)

#     println("="^70)
#     println("SANITY CHECK — shooting method / build_fAB")
#     println("L_box = $(round(L_box, digits=4)) r₀,  L/2 = $(round(L_box/2, digits=4)) r₀")
#     println("R0_test = $(round(R0_test, digits=4)) r₀")
#     println("="^70)

#     gr()

#     p_psi    = plot(xlabel = L"r/r_0", ylabel = L"\psi_b(r)\ /\ \psi_b(R_0)",
#                      title = "Dimer wavefunction (normalized at \$R_0\$)",
#                      legend = :topright, framestyle = :box, grid = true, gridalpha = 0.25)
#     p_uprime = plot(xlabel = L"r/r_0", ylabel = L"u'(r) = \psi_b'/\psi_b",
#                      title = "Logarithmic derivative", legend = :topright,
#                      framestyle = :box, grid = true, gridalpha = 0.25)
#     p_eb     = plot(xlabel = L"h/r_0", ylabel = L"\varepsilon_b", title = "Binding energy vs h",
#                      legend = false, framestyle = :box, grid = true, gridalpha = 0.25,
#                      marker = :circle)

#     eb_list = Float64[]

#     for h_test in h_vals_check
#         Δ = min(1e-4, h_test^(3/2) / 20.0)

#         result = build_fAB(h_test, R0_test, r_min, Δ, tol, nr0sq, N)

#         if result === nothing
#             println("h = $h_test  →  build_fAB returned `nothing` (node or ψ(R0)≈0)")
#             push!(eb_list, NaN)
#             continue
#         end

#         r_grid, psi, itp_u, itp_up, itp_upp, energy_b = result
#         push!(eb_list, energy_b)

#         u_prime_vals = [itp_up(r) for r in r_grid]

#         println("h = $(rpad(h_test, 5))  ε_b = $(round(energy_b, digits=6))  " *
#                 "ψ(r_min) = $(round(psi[1], digits=4))  ψ(R0) = $(round(psi[end], digits=4))  " *
#                 "min|ψ| = $(round(minimum(abs.(psi)), digits=8))")

#         plot!(p_psi, r_grid, psi, label = L"h = %$h_test", linewidth = 1.5)
#         plot!(p_uprime, r_grid, u_prime_vals, label = L"h = %$h_test", linewidth = 1.5)
#     end

#     plot!(p_eb, collect(h_vals_check), eb_list)

#     p_combined = plot(p_psi, p_uprime, p_eb;
#                        layout = (1, 3),
#                        size  = (1800, 550),
#                        left_margin  = 8Plots.mm,
#                        bottom_margin = 8Plots.mm)

#     plots_dir = normpath(joinpath(@__DIR__, "..", "data", "plots"))
#     mkpath(plots_dir)

#     display(p_combined)
#     savefig(p_combined, joinpath(plots_dir, "sanity_check_shooting_method.pdf"))

#     println("\n✓ Sanity-check plot saved to: $plots_dir")

#     # ════════════════════════════════════════════════════════════════
#     # SECOND CHECK — sweep R0 at FIXED h, replicating the exact range
#     # used by optimize_R0.jl: R0 ∈ [0.5h, L/2]. This tells you, before
#     # launching the full VMC sweep, which R0 points will come back
#     # `nothing` (node or ψ(R0)≈0) for a given h.
#     # ════════════════════════════════════════════════════════════════
#     h_fixed = 0.3
#     n_R0    = 20
#     R0_min  = 0.5 * h_fixed
#     R0_max  = L_box / 2
#     R0_vals_check = collect(LinRange(R0_min, R0_max, n_R0))

#     println("\n" * "="^70)
#     println("SECOND CHECK — R0 sweep at fixed h = $h_fixed")
#     println("R0 range : [$(round(R0_min,digits=4)), $(round(R0_max,digits=4))] r₀  ($(n_R0) points)")
#     println("="^70)

#     Δ_fixed = min(1e-4, h_fixed^(3/2) / 20.0)

#     p_psi_R0     = plot(xlabel = L"r/r_0", ylabel = L"\psi_b(r)\ /\ \psi_b(R_0)",
#                          title = "ψ_b(r) for several R0  (h = $h_fixed)",
#                          legend = false, framestyle = :box, grid = true, gridalpha = 0.25)
#     p_psiR0_raw  = plot(xlabel = L"R_0/r_0", ylabel = L"\psi_b(R_0)\ \mathrm{(unnormalized)}",
#                          title = "ψ_b(R0) before normalization", legend = false,
#                          framestyle = :box, grid = true, gridalpha = 0.25,
#                          marker = :circle, yscale = :log10)
#     p_uprime_R0  = plot(xlabel = L"R_0/r_0", ylabel = L"u'(R_0)",
#                          title = "u'(R0) — feeds kinetic energy", legend = false,
#                          framestyle = :box, grid = true, gridalpha = 0.25,
#                          marker = :circle)

#     R0_ok      = Float64[]
#     psiR0_raw  = Float64[]
#     uprime_R0  = Float64[]
#     R0_failed  = Float64[]

#     energy_b_fixed = find_energy_b(r_min, h_fixed, L_box/2, Δ_fixed, tol)
#     for R0 in R0_vals_check

#         y0 = initial_conditions(r_min, h_fixed, energy_b_fixed)
#         r  = r_min
#         n_steps_local = round(Int, (R0 - r_min) / Δ_fixed)
#         for _ in 1:n_steps_local
#             y0 = RK4(ode_rhs, r, y0, Δ_fixed, h_fixed, energy_b_fixed)
#             r += Δ_fixed
#         end
#         raw_psi_R0 = y0[1]

#         result = build_fAB(h_fixed, R0, r_min, Δ_fixed, tol, nr0sq, N)

#         if result === nothing
#             println("R0 = $(round(R0,digits=4))  →  nothing  (raw ψ(R0) = $(round(raw_psi_R0, sigdigits=4)))")
#             push!(R0_failed, R0)
#             continue
#         end

#         r_grid, psi, itp_u, itp_up, itp_upp, energy_b = result
#         push!(R0_ok, R0)
#         push!(psiR0_raw, abs(raw_psi_R0))
#         push!(uprime_R0, itp_up(r_grid[end]))

#         println("R0 = $(round(R0,digits=4))  raw ψ(R0) = $(round(raw_psi_R0, sigdigits=4))  " *
#                 "u'(R0) = $(round(itp_up(r_grid[end]), digits=4))")

#         plot!(p_psi_R0, r_grid, psi, label = L"R_0 = %$(round(R0,digits=2))", linewidth = 1.2, alpha = 0.7)
#     end

#     plot!(p_psiR0_raw, R0_ok, psiR0_raw)
#     plot!(p_uprime_R0, R0_ok, uprime_R0)

#     if !isempty(R0_failed)
#         println("\n⚠ $(length(R0_failed)) / $(n_R0) R0 points returned `nothing`: ",
#                 round.(R0_failed, digits=3))
#     else
#         println("\n✓ All $(n_R0) R0 points produced a valid f_AB")
#     end

#     p_combined2 = plot(p_psi_R0, p_psiR0_raw, p_uprime_R0;
#                         layout = (1, 3),
#                         size  = (1800, 550),
#                         left_margin  = 8Plots.mm,
#                         bottom_margin = 8Plots.mm)

#     display(p_combined2)
#     savefig(p_combined2, joinpath(plots_dir, "sanity_check_R0_sweep_h$(h_fixed).pdf"))

#     println("\n✓ R0-sweep sanity-check plot saved to: $plots_dir")
#     readline()
# end