module nDMolecularDynamicsSimulator
include("atoms.jl")

include("neighbor.jl")
include("potentials/atomic.jl")

include("input.jl")

function nDMD(inputfile::String)
    allSettings = loadfile(inputFile)
    run(allSettings)

end

end
