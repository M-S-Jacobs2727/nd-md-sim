export zeroForce!, LJCutForce!

function zeroForce!(atoms::Atoms, pairwiseParams::Array{<:Real})
    return
end

function LJCutForce!(atoms::Atoms, pairwiseParams::Array{<:Real})
    nLocal = atoms.nLocal
    f = zeros(Float64, size(atoms.forces, 1))
    for i=1:nLocal
        neighbors = atoms.neighbors[:, i]
        ri = atoms.positions[:, i]
        typeiParams = pairwiseParams[:, :, atoms.types[i]]
        for j in neighbors
            if j <= 0
                break
            end
            
            eps24, sigma6, rcutsq = typeiParams[:, atoms.types[j]]
            
            rj = atoms.positions[:, j]
            rdiff = rj - ri
            rsq = sum(rdiff.*rdiff)

            if rsq >= rcutsq
                continue
            end

            r6inv = 1/(rsq*rsq*rsq)
            f = rdiff * eps24/rsq*sigma6*r6inv * (2*sigma6*r6inv - 1)
            atoms.forces[:, i] .-= 0.5*f
            if j < nLocal + 1
                atoms.forces[:, j] .+= 0.5*f
            end
        end
    end
    return
end