export Langevin

struct Langevin <: Thermostat
    target::Function
    damping <: AbstractFloat
    function Langevin(target::Function, damping <: AbstractFloat)
        (damping > 0.0) || throw(ArgumentError("Argument damping must be positive."))
        applicable(target, 1.0) || throw(ArgumentError("Argument target should be a function taking one argument."))
        new(target, damping)
    end
    Langevin(target::AbstractFloat, damping <: AbstractFloat) = Langevin(damping, x -> target)
end

function postforce!(atoms::Atoms, thermostat::Langevin)

end

function integrate2!(atoms::Atoms, thermostat::Langevin)

end
