using Test

include("../src/nDMolecularDynamicsSimulator.jl")
using .nDMolecularDynamicsSimulator

@testset "LJ" begin
    # 2 atoms, 1 unit apart
    params = [24.0; 1.0; 6.25;;;]
    pos = [0 0; 0 0; 0 1]
    vel = [0 0; 0 0; 0 0]
    forces = [0 0; 0 0; 0 0]
    neighs = [2 1;]
    ids = [1, 2]
    types = [1, 1]
    masses = [1.0 1.0;]
    atoms = Atoms(2, 0, pos, vel, forces, neighs, ids, types, masses)

    LJCutForce!(atoms, params)
    @test isapprox(atoms.forces, [0 0; 0 0; -24.0 24.0])

    # 2 atoms, 2^(1/6) units apart (0 force)
    atoms.positions[3, 2] = 2^(1/6)
    atoms.forces[:] .= 0.0
    LJCutForce!(atoms, params)
    @test isapprox(atoms.forces, zeros(Float64, 3, 2))
    
    # 3 atoms, 2^(1/6) apart each
    atoms.nLocal = 3
    atoms.positions = [0 0 0; 0 0 0; -2^(1/6) 0 2^(1/6)]
    atoms.forces = zeros(Float64, 3, 3)
    atoms.types = [1, 1, 2]
    atoms.ids = [1, 2, 3]
    atoms.masses = ones(Float64, 1, 3)
    atoms.neighbors = [2 1 1; 3 3 2]
    params = zeros(Float64, 3, 2, 2)
    params[:, 1, 1] = [24.0, 1.0, 6.25]
    params[:, 1, 2] = [36.0, 1.0, 6.25]
    params[:, 2, 1] = [36.0, 1.0, 6.25]
    params[:, 2, 2] = [24.0, 1.0, 6.25]
    LJCutForce!(atoms, params)
    @test !isapprox(atoms.forces, [0 0 0; 0 0 0; 0 0 0])  # to complete later, if necessary        
end

@testset "Momentum" begin
    # 1 atom constant velocity
    params = zeros(Float64, 0, 0, 0)
    pos = zeros(Float64, 3, 1)
    vel = [0; 0; 1.0;;]
    forces = zeros(Float64, 3, 1)
    neighs = zeros(Int, 0, 0)
    ids = [1]
    types = [1]
    masses = [1.0;;]
    atoms = Atoms(1, 0, pos, vel, forces, neighs, ids, types, masses)
    settings = Settings(10, 0.1, 1, 1.0, 50, 50, params)
    subdomain = Subdomain(0, 1, nothing, zeros(3, 2), zeros(3, 2), zeros(Int, 3, 2))

    velocityVerlet!(atoms, settings, zeroForce!, subdomain)
    @test isapprox(atoms.positions, [0, 0, 1.0])

    # 2 atoms parallel
    settings = Settings(10, 0.1, 1, 50, 50, 50, [24.0; 1.0; 6.25;;;])
    atoms.nLocal = 2
    atoms.positions = [0 0; 0 0; 0 2^(1/6)]
    atoms.velocities = [1.0 1.0; 0 0; 0 0]
    atoms.forces = zeros(Float64, 3, 2)
    atoms.neighbors = [2 1;]
    atoms.ids = [1, 2]
    atoms.types = [1, 1]
    atoms.masses = [1.0 1.0;]
    velocityVerlet!(atoms, settings, LJCutForce!, subdomain)
    @test isapprox(atoms.positions, [1.0 1.0; 0 0; 0 2^(1/6)])
end