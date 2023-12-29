using LinearAlgebra: LinearAlgebra

export force!, maxneighbordist

abstract type AtomicPotential end

function force!(atoms::Atoms, potential::AtomicPotential, timestep::AbstractFloat) end

function maxneighbordist(potential::AtomicPotential, neighborskin::AbstractFloat) end

function symmetricify!(mat::Matrix{T <: float})
    (size(mat, 1) == size(mat, 2)) || throw(ArgumentError("Potential matrices must be square."))

    if istriu(mat)
        mat[:] .= Symmetric(mat, :U)
    elseif istril(mat)
        mat[:] .= Symmetric(mat, :L)
    elseif issymmetric(mat)
        throw(ArgumentError("Potential matrices must be either symmetric or triangular."))
    end
end
