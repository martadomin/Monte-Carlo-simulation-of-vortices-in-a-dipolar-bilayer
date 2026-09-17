# common/src/CommonCore.jl
#
# Include order follows the dependency chain — later files call functions
# defined in earlier ones:
#   utils               -> periodic geometry, RNG. No dependencies.
#   stats               -> pure statistics (blocking_statistics, detect_plateau).
#                          No dependencies.
#   jastrow_common      -> intra-layer/single-species Jastrow (u_AA family).
#                          No dependencies (besides the Bessels package).
#   jastrow_interlayer  -> inter-layer Jastrow (u_AB family), Stage2+3 only.
#                          No dependencies.
#   shooting_method     -> builds the interpolants jastrow_interlayer's u_AB
#                          consumes. No dependencies on the files above.
#   observables         -> uses get_periodic_difference (utils).
#   configuration       -> uses wrap_position, random_initial_config (utils).
#   trial_wavefunction  -> defines TrialWavefunction, the abstract type every
#                          logpsi/Δlogpsi/energy_estimators/init_observables/
#                          record_observables! method dispatches on.
#   dmc_kernel          -> diffusion_step uses wrap_position (utils).
#   metropolis          -> uses utils, configuration, and TrialWavefunction.
#   dmc_sampler         -> uses dmc_kernel, configuration, and TrialWavefunction.

include(joinpath(@__DIR__, "utils.jl"))
include(joinpath(@__DIR__, "stats.jl"))
include(joinpath(@__DIR__, "jastrow_common.jl"))
include(joinpath(@__DIR__, "jastrow_interlayer.jl"))
include(joinpath(@__DIR__, "shooting_method.jl"))
include(joinpath(@__DIR__, "observables.jl"))
include(joinpath(@__DIR__, "configuration.jl"))
include(joinpath(@__DIR__, "trial_wavefunction.jl"))
include(joinpath(@__DIR__, "dmc_kernel.jl"))
include(joinpath(@__DIR__, "metropolis.jl"))
include(joinpath(@__DIR__, "dmc_sampler.jl"))
