include("../src/JMD.jl")
import JMD

# set neighbor list settings
# set atomic potential
#   set potential parameters
# set constraints (constant volume, thermostat, barostat, etc.)
# read or create initial state
#   creation involves adding particles and velocities to regions of box
# run simulation

neighbor = JMD.Neighbor(0.3, 5)

potential = JMD.LJCut([1.0], [1.0], [2.5])

box = JMD.Box(JMD.BoxShape.ortho, [0.0 0.0 0.0; 10.0 10.0 10.0], fill(JMD.BoundaryCondition.periodic, 2, 3))
atoms = JMD.addatoms(box, density = 0.85, lattice = JMD.Lattice.bcc)
JMD.addvelocity!(atoms, temperature = 1.0)
tstat = JMD.Langevin(1.0, 10.0)
run = JMD.Run(0.005, 10_000)
sim = Simulation(atoms, neighbor, potential, box, [tstat])
JMD.verlet!(atoms, run, neighbor, potential, box, conditions = [tstat])

