# common/src/trial_wavefunction.jl

"""
    TrialWavefunction

Abstract supertype for all trial wavefunctions in this project. Every concrete
subtype represents one stage's physical system, and must implement:

    logpsi(trial::YourType, coords::NamedTuple, L) -> Float64
        Log of the trial wavefunction for the given configuration. `coords`
        is species-keyed, e.g. (A=(x=[...],y=[...]), B=(x=[...],y=[...])).

    Δlogpsi(trial::YourType, coords::NamedTuple, species::Symbol, id::Int,
            x_new::Float64, y_new::Float64, L) -> Float64
        Change in logpsi from moving particle `id` of `species` in `coords`
        to (x_new, y_new). Only needs to be correct for a single-particle
        move (this is what move_one_part/tune_delta use).

    energy_estimators(trial::YourType, coords::NamedTuple, L) -> NamedTuple
        Fields: drift (NamedTuple, coords-shaped — e.g. drift.A.x, drift.B.x
        — one entry per particle per species, matching coords' own keys),
        E_total_local, E_total_drift, E_total_laplacian (three independent
        estimators of the same total energy, kinetic+interaction[+external];
        should agree within statistical error — a persistent mismatch signals
        a bug), E_kinetic (kinetic energy alone, local/std estimator),
        E_interaction (pairwise interaction energy alone — any one-body
        external potential a stage has, e.g. a vortex phase term, is folded
        into E_total_* but NOT into this field).

    init_observables(trial::YourType, num_bins::Int) -> Any
        Builds whatever observable-accumulator object this stage tracks.
        Shape is entirely this stage's choice — not part of the contract
        beyond "something record_observables! can update."

    record_observables!(trial::YourType, obs, coords::NamedTuple, L)
        Accumulates one configuration's contribution into `obs`.

metropolis.jl and dmc_sampler.jl call all five generically, without knowing
which concrete type they're handed — adding a new physical system means
adding a new subtype and these five methods here, not touching the samplers.
"""
abstract type TrialWavefunction end