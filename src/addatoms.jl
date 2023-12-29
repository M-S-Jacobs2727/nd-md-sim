export addatoms, Lattice

@enum Lattice begin
    random
    fcc
    bcc
    hcp
    sq
    sq2
    cu
end

function addatoms(box::Box; density <: float, lattice::Lattice = Lattice.random)::Atoms

end

function addatoms(box::Box; spacing <: float, lattice::Lattice = Lattice.cu)::Atoms
    (lattice == Lattice.random) && throw(ArgumentError(""))
end
