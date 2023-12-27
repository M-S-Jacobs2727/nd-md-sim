export Atoms

struct Atoms{T1 <: Unsigned, T2 <: AbstractFloat}
    natoms::T1
    dim::T1
    ids::Vector{T1}
    types::Vector{T1}
    positions::Matrix{T2}
    velocities::Matrix{T2}
    forces::Matrix{T2}
    neighbors::Vector{Vector{T1}}
end

Atoms{IntType <: Unsigned, FloatType <: AbstractFloat}(natoms::IntType, dim::IntType = 3) = Atoms{IntType, FloatType}(
    natoms,
    dim,
    zeros(IntType, natoms),
    zeros(IntType, natoms),
    zeros(FloatType, dim, natoms),
    zeros(FloatType, dim, natoms),
    zeros(FloatType, dim, natoms),
    fill(Vector{IntType}(), natoms),
)
