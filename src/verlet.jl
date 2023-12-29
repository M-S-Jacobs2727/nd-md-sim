function halfvelocity!{T}(atoms::Atoms, timestep::T) where T <: AbstractFloat
    @. atoms.velocities += (0.5 * timestep) * atoms.forces / atoms.masses
end

function positionstep!{T}(atoms::Atoms, timestep::T) where T <: AbstractFloat
    @. atoms.positions += timestep * (atoms.velocities + (0.5 * timestep) * atoms.forces / atoms.masses)
end

function verlet!(sim::Simulation; timestep::AbstractFloat, steps::Integer)
    (timestep <= 0) && throw(DomainError("Timestep must be positive."))
    (steps <= 0) && throw(DomainError("Number of steps must be positive."))

    maxdist = maxneighbordist(potential, sim.neighbor.cutoff)

    pbc!(sim.atoms, sim.box)
    reneighbor!(sim.atoms, maxdist)

    sim.atoms.forces .= 0.0
    force!(sim.atoms, sim.potential)

    for i ∈ 1:steps
        for condition in sim.conditions
            integrate1!(sim.atoms, condition)
        end

        halfvelocity!(sim.atoms, timestep)
        positionstep!(sim.atoms, timestep)

        if (i - 1) % sim.neighbor.updateevery == 0 && needs_reneighbor(sim.atoms, maxdist)
            pbc!(sim.atoms, sim.box)
            reneighbor!(sim.atoms, maxdist)
        end

        sim.atoms.forces .= 0.0
        force!(sim.atoms, sim.potential)

        halfvelocity!(sim.atoms, timestep)
    end
end
