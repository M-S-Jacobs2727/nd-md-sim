export Neighbor

struct Neighbor
    skin::AbstractFloat
    updateevery::Unsigned
end

function needs_reneighbor(atoms::Atoms, maxdist::AbstractFloat)::Bool
end

function reneighbor!(atoms::Atoms, maxdist::AbstractFloat)
    # for each atom in writeview
    #   for each atom in atom's neighborset
    #       if distance > maxdist: remove from neighborset
    #   for each atom in localview and ghostview
    #       if not in neighborset and distance <= maxdist: add to neighborset
end
