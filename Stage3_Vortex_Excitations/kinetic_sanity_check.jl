"""
test_kinetic_identity.jl

Standalone diagnostic: checks the identity
    K_i = (1/2) F_i^2 = -(1/4) T_i
for a SINGLE pair of particles, comparing:
  (a) the analytic u', u'' from jastrow.jl / shooting_method.jl
  (b) finite-difference derivatives of u and u' computed directly
      from u_AA(r,...) and u_AB(r,...) (i.e. from f, not from f', f'')

If (a) and (b) disagree for AA  -> bug is in u_AA_prime/u_AA_second (old, Stage I code)
If (a) and (b) disagree for AB  -> bug is in shooting_method.jl's u_prime/u_doubleprime
If (a) and (b) AGREE but the energy.jl combination still doesn't satisfy the
   K = 1/2 F^2 = -1/4 T identity in a real multi-particle config, the bug is
   in how energy_estimators() assembles the per-particle sums.

Run directly: julia .\\src\\test_kinetic_identity.jl
(adjust the include paths below to your actual src/ layout)
"""

include(normpath(joinpath(@__DIR__, "src",  "jastrow.jl")))
include(normpath(joinpath(@__DIR__, "src", "shooting_method.jl")))

eps = 1e-6

println("="^70)
println("TEST 1 — u_AA, u_AA_prime, u_AA_second: analytic vs finite difference")
println("="^70)

L = 7.745966692414834
R_match = 0.8349980438651612
Constants = calculate_constants(L, R_match)

global max_err_prime_AA  = 0.0
global max_err_second_AA = 0.0

for r in [0.1, 0.3, 0.5, 0.835, 1.0, 1.5, 2.0, 3.0, 3.5]
    u0 = u_AA(r, R_match, L, Constants)
    up = u_AA(r + eps, R_match, L, Constants)
    um = u_AA(r - eps, R_match, L, Constants)

    fd_prime  = (up - um) / (2eps)
    fd_second = (up - 2u0 + um) / eps^2

    an_prime  = u_AA_prime(r, R_match, L, Constants)
    an_second = u_AA_second(r, R_match, L, Constants)

    err_p = abs(fd_prime - an_prime)
    err_s = abs(fd_second - an_second)
    global max_err_prime_AA  = max(max_err_prime_AA, err_p)
    global max_err_second_AA = max(max_err_second_AA, err_s)

    println("r=$(rpad(r,5))  u'_an=$(round(an_prime,digits=6))  u'_FD=$(round(fd_prime,digits=6))  " *
            "Δ=$(round(err_p,sigdigits=3))  ||  u''_an=$(round(an_second,digits=4))  " *
            "u''_FD=$(round(fd_second,digits=4))  Δ=$(round(err_s,sigdigits=3))")
end

println("\nMax error AA:  u' -> $max_err_prime_AA   u'' -> $max_err_second_AA")
println(max_err_prime_AA < 1e-3 && max_err_second_AA < 1e-1 ?
        "✓ AA derivatives look consistent" : "✗ AA derivatives DISAGREE with finite differences")

println("\n" * "="^70)
println("TEST 2 — u_AB, u_AB_prime, u_AB_second: analytic vs finite difference")
println("="^70)

r_min   = 1e-6
tol     = 1e-10
nr0sq   = 1.0
N       = 60
h_test  = 0.78           # match the case you plotted
L_box   = sqrt(N / nr0sq)
R0_test = 0.3             # near your reported optimum for h=0.78

Δ = min(1e-4, h_test^(3/2) / 20.0)
result = build_fAB(h_test, R0_test, r_min, Δ, tol, nr0sq, N)

max_err_prime_AB  = 0.0
max_err_second_AB = 0.0

if result === nothing
    println("build_fAB returned nothing for h=$h_test, R0=$R0_test — adjust R0_test and rerun.")
else
    let max_err_prime_AB = 0.0, max_err_second_AB = 0.0
        r_grid, psi, u_prime_grid, u_doubleprime_grid, energy_b = result
        println("h=$h_test  R0=$R0_test  ε_b=$(round(energy_b,digits=6))")

        # sample a few interior points of the grid. For u_AB' the finite-difference
        # check is meaningful because u_AB itself is a smooth interpolant.
        # For u_AB_second, direct finite differences of u_AB are misleading because
        # u_AB is piecewise linear in log(psi) between grid nodes. Instead, compare
        # the stored u_doubleprime grid against a finite difference of u_prime on
        # the same ODE grid.
        test_idxs = round.(Int, LinRange(5, length(r_grid)-5, 8))

        for idx in test_idxs
            r = r_grid[idx]

            u0 = u_AB(r, R0_test, r_grid, psi)
            up = u_AB(r + eps, R0_test, r_grid, psi)
            um = u_AB(r - eps, R0_test, r_grid, psi)

            fd_prime = (up - um) / (2eps)

            an_prime  = u_AB_prime(r, R0_test, r_grid, u_prime_grid)
            an_second = u_AB_second(r, R0_test, r_grid, u_doubleprime_grid)

            err_p = abs(fd_prime - an_prime)
            max_err_prime_AB  = max(max_err_prime_AB, err_p)

            fd_second_grid = (u_prime_grid[idx+1] - u_prime_grid[idx-1]) /
                     (r_grid[idx+1] - r_grid[idx-1])
            err_s = abs(fd_second_grid - u_doubleprime_grid[idx])
            max_err_second_AB = max(max_err_second_AB, err_s)

            println("r=$(round(r,digits=4))  u'_an=$(round(an_prime,digits=6))  " *
                "u'_FD=$(round(fd_prime,digits=6))  Δ=$(round(err_p,sigdigits=3))  ||  " *
                "u''_grid=$(round(u_doubleprime_grid[idx],digits=4))  " *
                "u''_FD(grid)=$(round(fd_second_grid,digits=4))  " *
                "Δ=$(round(err_s,sigdigits=3))")
        end

        println("\nMax error AB:  u' -> $max_err_prime_AB   u'' -> $max_err_second_AB")
        println(max_err_prime_AB < 1e-2 && max_err_second_AB < 1e-2 ?
            "✓ AB derivatives look consistent" : "✗ AB derivatives need closer inspection")
    end
end

println("\n" * "="^70)
println("TEST 3 — single-pair decomposition and on-shell Schrödinger residual")
println("="^70)
println("For a single AB pair at separation r, with no other particles,")
println("check consistency of T + F^2 with (V_AB - ε_b) from the bound-state ODE.")
println("="^70)

if result !== nothing
    r_grid, psi, u_prime_grid, u_doubleprime_grid, energy_b = result
    for r_test in [0.05, 0.1, 0.15, 0.2, 0.25]
        up = u_AB_prime(r_test, R0_test, r_grid, u_prime_grid)
        u2 = u_AB_second(r_test, R0_test, r_grid, u_doubleprime_grid)

        F_sq = up^2          # |∇u|^2 in radial form (single pair, 1D radial gradient squared)
        T    = u2 + up/r_test   # full 2D radial Laplacian of u
        V    = V_AB(r_test, h_test)
        rhs  = V - energy_b
        schr_residual = (T + F_sq) - rhs

        K_full  = -0.5 * (F_sq + T)
        K_drift =  0.5 * F_sq
        K_lap   = -0.25 * T

        println("r=$r_test  K_full=$(round(K_full,digits=6))  " *
                "K_drift=$(round(K_drift,digits=6))  K_lap=$(round(K_lap,digits=6))  " *
                "T+F²=$(round(T+F_sq,digits=6))  V-ε=$(round(rhs,digits=6))  " *
                "res=$(round(schr_residual,sigdigits=3))")
    end
end