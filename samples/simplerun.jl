include("../src/nDMolecularDynamicsSimulator.jl")
import nDMolecularDynamicsSimulator as ndmds

# set neighbor list settings
# set atomic potential
#   set potential parameters
# set constraints (constant volume, thermostat, barostat, etc.)
# read or create initial state
#   creation involves adding particles and velocities to regions of box
# run simulation

neighbor = ndmds.Neighbor(0.3, 5)

potential = ndmds.LJCut([1.0], [1.0], [2.5])

box = ndmds.Box(ndmds.BoxShape.ortho, [0.0 0.0 0.0; 10.0 10.0 10.0], fill(ndmds.BoundaryCondition.periodic, 2, 3))
atoms = ndmds.addatoms(box, density = 0.85, lattice = ndmds.Lattice.bcc)
ndmds.addvelocity!(atoms, temperature = 1.0)
tstat = ndmds.Langevin(1.0, 10.0)
run = ndmds.Run(0.005, 10_000)
sim = Simulation(atoms, neighbor, potential, box, [tstat])
ndmds.verlet!(atoms, run, neighbor, potential, box, conditions = [tstat])

