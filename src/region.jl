using LinearAlgebra
import Base: in

abstract type Manifold end # (n-1)-dimensional manifold in n-dimensional space

struct Plane{T <: AbstractFloat} <: Manifold
    refpoint::Vector{T}
    normal::Vector{T}
    function Plane(refpoint::Vector{T}, normal::Vector{T}) where T <: AbstractFloat
        (length(refpoint) == length(normal)) ||
            throw(DimensionMismatch("refpoint and normal must be the same dimension (length)."))
        normalize!(normal)
        new(refpoint, normal)
    end
end

struct Sphere{T <: AbstractFloat} <: Manifold
    center::Vector{T}
    radius::T
end

struct Subspace{T <: AbstractFloat} # m-dimensional manifold in n-dimensional space, m ∈ 1:(n-1)
    refpoint::Vector{T}
    basis::Matrix{T}
    function Subspace(refpoint::Vector{T}, basis::Matrix{T}) where T <: AbstractFloat
        (size(refpoint, 1) == size(basis, 1)) ||
            throw(DimensionMismatch("refpoint and basis matrix must be same size in first dimension"))

        (size(basis, 1) > size(basis, 2)) ||
            throw(DimensionMismatch("basis must form a proper subspace (more rows than columns)"))

        r = rank(basis)
        (r == size(basis, 2)) ||
            throw(ArgumentError("basis must consist of linearly independent columns."))

        for c ∈ eachcol(basis)
            normalize!(c)
        end
        new(refpoint, basis)
    end
end

struct Cylinder{T <: AbstractFloat} <: Manifold
    axis::Subspace{T}
    radius::T
end

function in(point::Vector{T}, plane::Plane{T}) where T <: AbstractFloat
    return dot(point .- plane.refpoint, plane.normal) > zero(T)
end

function in(point::Vector{T}, sphere::Sphere{T}) where T <: AbstractFloat
    return norm(point .- sphere.center) < sphere.radius
end

function in(point::Vector{T}, cyl::Cylinder{T}) where T <: AbstractFloat
    r = point .- cyl.subspace.refpoint
    for b ∈ eachcol(cyl.subspace.basis)
        r .-= b .* dot(b, r)
    end
    return norm(r) < cyl.radius
end

struct Region
    manifolds::Set{Manifold} = Set(Manifold[])
end

function Region(box::Box{T}) where T <: AbstractFloat
    manifolds = Set(Manifold[])

    refpoint = zeros(T, size(box.bounds, 1))
    norm = zeros(T, size(box.bounds, 1))

    for (i, (lo, hi)) ∈ enumerate(eachrow(box.bounds))
        refpoint[i] = lo
        norm[i] = one(T)
        union!(manifolds, Plane{T}(refpoint, norm))

        refpoint[i] = hi
        norm[i] *= -1
        union!(manifolds, Plane{T}(refpoint, norm))

        refpoint[i] = zero(T)
        norm[i] = zero(T)
    end

    Region(manifolds)
end

function in(point::Vector{T}, region::Region{T})
    for m ∈ region.manifolds
        (point ∉ m) && return false
    end
    return true
end
