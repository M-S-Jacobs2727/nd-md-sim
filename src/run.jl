export run!, RunSettings

struct RunSettings
    timestep::float
    steps::Integer
end

function run!(atoms::Atoms, run::RunSettings, neighbor::NeighborSettings)
end