struct Atoms{U <: Unsigned, F <: AbstractFloat}
    natoms::U
    dim::U
    ids::Vector{U}
    types::Vector{U}
    masses::Vector{F}
    positions::Matrix{F}
    velocities::Matrix{F}
    forces::Matrix{F}
    neighbors::Vector{Set{U}}
    function Atoms{U, F}(
        natoms::U,
        dim::U,
        ids::Vector{U},
        types::Vector{U},
        masses::Vector{F},
        positions::Matrix{F},
        velocities::Matrix{F},
        forces::Matrix{F},
        neighbors::Vector{Set{U}},
    ) where U <: Unsigned where F <: AbstractFloat
        (1 <= dim <= 6) || throw(ArgumentError("Number of dimensions must be between 1 and 6 inclusive."))

        (length(ids) == natoms) && (length(types) == natoms) && (length(masses) == natoms) && (length(neighbors) == natoms) ||
            throw(ArgumentError("Size of ids, types, masses and neighbors must be equal to natoms."))

        (size(positions) == (dim, natoms)) && (size(velocities) == (dim, natoms)) && (size(forces) == (dim, natoms)) ||
            throw(ArgumentError("Size of positions, velocities, and forces arrays must be equal to dim x natoms."))

        new(natoms, dim, ids, types, masses, positions, velocities, forces, neighbors)
    end
end

Atoms{U, F}(natoms::Integer, dim::Integer = 3) where U <: Unsigned where F <: AbstractFloat =
    Atoms{U, F}(
        U(natoms),
        U(dim),
        zeros(U, U(natoms)),
        zeros(U, U(natoms)),
        ones(F, U(natoms)),
        zeros(F, U(dim), U(natoms)),
        zeros(F, U(dim), U(natoms)),
        zeros(F, U(dim), U(natoms)),
        fill(Set{U}(), U(natoms)),
    )
