module nDMolecularDynamicsSimulator
export Atoms, Settings, main, velocityVerlet!
import MPI
import Combinatorics: partitions

include("subdomain.jl")

mutable struct Atoms  # per-proc
    nLocal::Int
    nGhost::Int
    positions::Matrix{Float64}
    velocities::Matrix{Float64}
    forces::Matrix{Float64}
    neighbors::Matrix{Int}
    ids::Vector{Int}
    types::Vector{Int}
    masses::Matrix{Float64}
    images::Matrix{Int}
end
include("pairwisePotentials.jl")
struct Settings  # global
    timesteps::Int
    dt::Float64
    nAtomTypes::Int
    neighborCutoff::Float64
    neighborUpdateEvery::Int
    maxNeighborsPerAtom::Int
    pairwiseParams<:pairwisePotential
end
include("data.jl")
include("input.jl")

function updateNeighbors!(atoms::Atoms, neighborCutoff::Float64, subdomain::Subdomain)
end

function velocityVerlet!(atoms::Atoms, settings::Settings, subdomain::Subdomain)
    timesteps = settings.timesteps
    dt = settings.dt
    neighborCutoff = settings.neighborCutoff
    neighborUpdateEvery = settings.neighborUpdateEvery
    pairwiseParams = settings.pairwiseParams
    
    atoms.forces[:] .= 0.0
    updateNeighbors!(atoms, neighborCutoff, subdomain)
    pairwiseForce!(atoms, pairwiseParams)
    for i=1:timesteps
        atoms.positions += atoms.velocities*dt + 0.5*dt*dt*atoms.forces ./ atoms.masses
        atoms.velocities += 0.5*dt*atoms.forces ./ atoms.masses
        if (i-1)%neighborUpdateEvery == 0
            updateNeighbors!(atoms, neighborCutoff, subdomain)
        end
        atoms.forces[:] .= 0.0
        pairwiseForce!(atoms, pairwiseParams)
        atoms.velocities += 0.5*dt*atoms.forces ./ atoms.masses
    end
end

# function balanceProcs!(subdomain::Subdomain, settings::Settings, atoms::Atoms)
#     localPositions = atoms.positions[:, 1:atoms.nLocal]
#     nAtoms = 0
#     nAtoms = MPI.Reduce(atoms.nLocal, MPI.SUM, 0, subdomain.comm)
#     q, r = divrem(nAtoms, subdomain.nProcs)

#     for (dim, gridLength) in enumerate(settings.nProcsPerDim)

#     end

# end

function main(datafile::String, inputfile::String, nDim::Int)
    MPI.Init()
    comm = MPI.COMM_WORLD
    me = MPI.Comm_rank(comm)
    nProcs = MPI.Comm_size(comm)
    if me == 0
        atoms::Atoms, boxBounds = readDatafile(datafile, nDim)
    else
        atoms = Atoms(nDim)
        boxBounds = zeros(Float64, nDim, 2)
    end
    atoms = MPI.bcast(atoms, 0, comm)  # don't duplicate all the atoms, distribute them in-place
    MPI.Bcast!(boxBounds, 0, comm)
    subdomain::Subdomain = setupSubdomain(me, nProcs, comm, boxBounds)
    distributeAtoms!(atoms, subdomain)
    parseInput(inputfile, atoms, subdomain)
end

end # module nDMolecularDynamicsSimulator
