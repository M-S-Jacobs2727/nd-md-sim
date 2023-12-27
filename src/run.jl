export run!, Run

struct Run
    timestep <: float
    steps <: Unsigned
    function Run(timestep <: float, steps <: Unsigned)
        (timestep <= 0) && throw(DomainError("Timestep must be positive."))
        new(timestep, steps)
    end
end

function run!(atoms::Atoms, run::Run, neighbor::NeighborSettings)
end