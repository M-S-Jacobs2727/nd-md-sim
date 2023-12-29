struct Simulation
    atoms::Atoms
    neighbor::Neighbor
    potential::AtomicPotential
    box::Box
    dimensions::Unsigned = 3
    conditions::Vector{SimCondition} = SimCondition[]
    Simulation(
        atoms::Atoms,
        neighbor::Neighbor,
        potential::AtomicPotential,
        box::Box,
        dimensions::Unsigned = 3,
        conditions::SimCondition...,
    ) = Simulation(atoms, neighbor, potential, box, dimensions, conditions)
end
