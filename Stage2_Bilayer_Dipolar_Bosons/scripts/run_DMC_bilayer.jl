include(normpath(joinpath(@__DIR__, "..", "src", "jastrow.jl")))
include(normpath(joinpath(@__DIR__, "..", "src", "shooting_method.jl")))
include(normpath(joinpath(@__DIR__, "..", "src", "energy.jl")))
include(normpath(joinpath(@__DIR__, "..", "src", "utils.jl")))
include(normpath(joinpath(@__DIR__, "..", "src", "observables.jl")))
include(normpath(joinpath(@__DIR__, "..", "src", "dmc.jl")))
include(normpath(joinpath(@__DIR__, "..", "src", "metropolis.jl")))


# Run single DMC at h = 1.0
N = 60
nr0sq = 1.0
h_dmc   = 0.3
N_half  = N ÷ 2


R0_opt      = 0.7281208690869945
energy_b    = -11.661903206745434
E_opt       = -0.6262732051697336 * 60
tail_opt    = 0.807532720541924

r_min     = 1e-6
Δ_shoot   = 1e-4
tol_shoot = 1e-10

L         = sqrt(N / nr0sq)
R_match   = 0.8349980438651612   # from Stage I
Constants = calculate_constants(L, R_match)

# Build fAB at optimal R0 for this h
result = build_fAB(h_dmc, R0_opt, r_min, Δ_shoot, tol_shoot, nr0sq, N)
r_grid_dmc, psi_dmc, itp_u_dmc, itp_up_dmc, itp_upp_dmc, energy_b_dmc = result

# Initial config from VMC
x_coord, y_coord = random_initial_config(N, L, "Uniform")
xA_init = x_coord[1:N_half];     yA_init = y_coord[1:N_half]
xB_init = x_coord[N_half+1:end]; yB_init = y_coord[N_half+1:end]

# Tune delta first
delta_dmc, xA_init, yA_init, xB_init, yB_init = tune_delta(
    xA_init, yA_init, xB_init, yB_init,
    L, R_match, Constants, R0_opt, itp_u_dmc;
    target_ratio = 0.5, num_tune_steps = 5000)

# DMC parameters
Δτ = 0.0004 * 0.1       # from δ = sqrt(2D Δτ), D=0.5
num_walkers = 400
num_target  = 400
num_steps   = 10^5
E_ref_init  = E_opt * N           # use VMC energy as starting reference

E_dmc, E_dmc_err, E_history = dmc(
    xA_init, yA_init, xB_init, yB_init,
    num_walkers, N, num_steps, Δτ,
    L, h_dmc, R_match, Constants, R0_opt,
    itp_up_dmc, itp_upp_dmc,
    E_ref_init, num_target;
    num_equil    = num_steps ÷ 5,
    quadratic    = true,
    plot_energy  = true
)

println("DMC E/N = $(round(E_dmc/N, digits=6)) ± $(round(E_dmc_err/N, digits=6))")