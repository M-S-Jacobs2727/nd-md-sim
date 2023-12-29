export LJCut

struct LJCut{T <: AbstractFloat} <: AtomicPotential
    sigma::Matrix{T}
    epsilon::Matrix{T}
    rcut::Matrix{T}
    sigma6::Matrix{T}
    epsilon12::Matrix{T}
    rcutsq::Matrix{T}
    function LJCut{T}(
        sigma::Matrix{T},
        epsilon::Matrix{T},
        rcut::Matrix{T},
        sigma6::Matrix{T},
        epsilon12::Matrix{T},
        rcutsq::Matrix{T},
    ) where T <: AbstractFloat
        (size(sigma) == size(epsilon) == size(rcut) == size(sigma6) == size(epsilon12) == size(rcutsq)) ||
            throw(ArgumentError("Potential matrices must be the same sizes."))

        for mat in [sigma, epsilon, rcut, sigma6, epsilon12, rcutsq]
            symmetricify!(mat)
        end
        new(sigma, epsilon, rcut, sigma6, epsilon12, rcutsq)
    end
end

LJCut(sigma, epsilon, rcut) = LJCut(sigma, epsilon, rcut, sigma .^ 6, 12 .* epsilon, rcut .^ 2)

function force!(atoms::Atoms, params::LJCut)
    for i in 1:size(atoms.neighbors, 2), j in atoms.neighbors[i]
        rdiff::Vector{float} = atoms.positions[:, j] - atoms.positions[:, i]
        rsq::float = sum(rdiff .^ 2)

        rcutsq = params.rcutsq[atoms.types[i], atoms.types[j]]
        (rsq > rcutsq) && continue

        epsilon12 = params.epsilon12[atoms.types[i], atoms.types[j]]
        sigma6 = params.sigma6[atoms.types[i], atoms.types[j]]

        r6inv = 1 / rsq^3

        # this is half the force; the full force would be 24*epsilon instead of 12*epsilon
        atoms.forces[:, i] .-= rdiff / rsq * epsilon12 * sigma6 * r6inv * (2 * sigma6 * r6inv - 1)
    end
end

function maxneighbordist{T <: AbstractFloat}(potential::LJCut{T}, neighborskin::{T})
    potential.rcut + neighborskin
end
