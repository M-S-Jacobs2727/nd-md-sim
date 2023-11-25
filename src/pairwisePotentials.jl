export pairwiseForce!, pairwisePotential, LJCutPotential, zeroPotential, nmCutPotential

abstract type pairwisePotential end

struct zeroPotential<:pairwisePotential end

struct LJCutPotential<:pairwisePotential
    nParams::Int
    epsilon24::Matrix{Float64}
    sigma6::Matrix{Float64}
    rcut2::Matrix{Float64}
    function LJCutPotential(;epsilon::Matrix{Float64},
            sigma::Matrix{Float64},
            rcut::Matrix{Float64})
        sigma3 = sigma*sigma*sigma
        new(24*epsilon, sigma3*sigma3, rcut*rcut)
    end
end

struct nmCutPotential<:pairwisePotential
    nParams::Int
    epsilon::Matrix{Float64}
    sigma::Matrix{Float64}
    rcut::Matrix{Float64}
    n::Int
    m::Int
end

function pairwiseForce!(atoms::Atoms, params::zeroPotential)
    return
end

function pairwiseForce!(atoms::Atoms, params::LJCutPotential)
    nLocal = atoms.nLocal
    f = zeros(Float64, size(atoms.forces, 1))
    for i=1:nLocal
        neighbors = atoms.neighbors[:, i]
        ri = atoms.positions[:, i]
        for j in neighbors
            if j <= 0
                break
            end
            
            epsilon24 = params.epsilon24[atoms.types[j], atoms.types[i]]
            sigma6 = params.sigma6[atoms.types[j], atoms.types[i]]
            rcut2 = params.rcut2[atoms.types[j], atoms.types[i]]
            
            rj = atoms.positions[:, j]
            rdiff = rj - ri
            rsq = sum(rdiff.*rdiff)

            if rsq >= rcut2
                continue
            end

            r6inv = 1/(rsq*rsq*rsq)
            f = rdiff * epsilon24/rsq*sigma6*r6inv * (2*sigma6*r6inv - 1)
            atoms.forces[:, i] .-= 0.5*f
            if j < nLocal + 1
                atoms.forces[:, j] .+= 0.5*f
            end
        end
    end
    return
end