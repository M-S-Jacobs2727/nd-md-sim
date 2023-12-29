export LJCut

struct LJCut <: AtomicPotential
    sigma::float
    epsilon::float
    rcut::float
    sigma6::float
    epsilon24::float
end

LJCut(sigma, epsilon, rcut) = LJCut(sigma, epsilon, rcut, sigma^6, 24 * epsilon)

function force!(atoms::Atoms, params::LJCut)
    for i in 1:size(atoms.neighbors, 2), j in atoms.neighbors[i]
        rdiff::Vector{float} = atoms.positions[:, j] - atoms.positions[:, i]
        rsq::float = sum(rdiff .^ 2)
        r6inv = 1 / rsq^3
        atoms.forces[:, i] .-= 0.5 * rdiff * params.epsilon24 / rsq * params.sigma6 * r6inv * (2 * params.sigma6 * r6inv - 1)
    end
end
