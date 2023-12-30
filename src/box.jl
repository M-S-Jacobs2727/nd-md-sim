export Box, BoxShape, BoundaryCondition, pbc!

@enum BoxShape begin
    ortho
    triclinic
end

@enum BoundaryCondition begin
    periodic
    fixed
    minimum
    shrink
end

struct Box{T1 <: AbstractFloat}
    shape::BoxShape
    bounds::Matrix{T1}
    pbc::Matrix{BoundaryCondition}
    function Box(shape::BoxShape, bounds::Matrix{T}, pbc::Matrix{BoundaryCondition}) where T <: AbstractFloat
        any(bounds[:, 1] .>= bounds[:, 2]) && throw(
            ArgumentError("Box bounds must be a 2xn matrix with bounds[:, 1] < bounds[:, 2]'."),
        )
        (size(bounds) == size(pbc)) || throw(ArgumentError("Sizes of bounds and pbc must match."))
        (size(bounds, 2) == 2) || throw(ArgumentError("Second dimension of bounds and pbc must be 2."))
        for dim ∈ eachrow(pbc)
            # If exactly one of the pbcs in a given dimension is periodic, throw
            ((dim[1] == periodic) ⊻ (dim[2] == periodic)) && throw(
                ArgumentError("Box dimensions with periodic boundaries must be periodic on both sides."),
            )
        end

        new{T}(shape, bounds, pbc)
    end
end

# function Box{T <: AbstractFloat}(shape::AbstractString, bounds::Vector{T}, pbc::AbstractString)
#     if shape == "ortho"
#         _shape = BoxShape.ortho
#     elseif shape == "triclinic"
#         _shape = BoxShape.triclinic
#     else
#         throw(ArgumentError("Unknown box shape; must be either 'ortho' or 'triclinic'."))
#     end

#     _bounds = reshape(bounds, 2, :)

#     _pbc = fill(BoundaryCondition.periodic, 2, dim)
#     for i ∈ 1:2dim
#         if pbc[i] == 'f'
#             _pbc[i] = BoundaryCondition.fixed
#         elseif pbc[i] == 'm'
#             _pbc[i] = BoundaryCondition.minimum
#         elseif pbc[i] == 's'
#             _pbc[i] = BoundaryCondition.shrink
#         elseif pbc[i] != 'p'
#             throw(ArgumentError("Box pbc must be defined with the characters p, m, f, or s."))
#         end
#     end

#     Box{T}(_shape, _bounds, _pbc)
# end
