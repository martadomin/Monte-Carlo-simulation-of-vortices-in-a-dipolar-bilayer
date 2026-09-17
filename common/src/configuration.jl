"""
    copy_coords(coords) -> NamedTuple

Deep copy of a coords NamedTuple — new arrays, not just a new NamedTuple
wrapping the same ones. Used at every point external `coords` enters a
function that's going to mutate its working copy.
"""
copy_coords(coords::NamedTuple) = map(c -> (x=copy(c.x), y=copy(c.y)), coords)
"""
    init_random_config(n_species, L; distribution="Uniform") -> NamedTuple

Builds a fresh coords NamedTuple with `n_species[s]` particles of species
`s`, uniformly placed in the box, for every species named in `n_species`.
`n_species` is (A=30, B=30)-shaped (species => particle count, plain Ints);
the return value is (A=(x=...,y=...), B=(x=...,y=...))-shaped — this is
the one function that turns the first shape into the second.
"""
function init_random_config(n_species::NamedTuple, L::Float64; distribution::AbstractString="Uniform")
    return map(n_species) do n
        x, y = random_initial_config(n, L, distribution)
        (x=x, y=y)
    end
end

species_names(coords::NamedTuple) = keys(coords)
num_particles(coords::NamedTuple) = sum(length(c.x) for c in coords)

"""
    random_particle(coords) -> (species::Symbol, local_id::Int)

Picks one particle uniformly at random across all species.
"""
function random_particle(coords::NamedTuple)
    names  = species_names(coords)
    counts = [length(coords[s].x) for s in names]
    idx    = rand(1:sum(counts))
    cum    = 0
    for (s, c) in zip(names, counts)
        idx <= cum + c && return s, idx - cum
        cum += c
    end
end