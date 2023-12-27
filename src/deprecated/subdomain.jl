export Subdomain, setupSubdomain

mutable struct Subdomain  # per-proc
    me::Int  # MPI rank, 0-based
    nProcs::Int  # number of MPI ranks
    comm::MPI.Comm  # MPI communicator, usually COMM_WORLD
    localBounds::Matrix{Rational}  # subdomain bounds in scaled [0,1] coordinates
    boxBounds::Matrix{Float64}  # whole box bounds, standard coordinates
    neighborProcs::Matrix{Int}  # MPI rank of the 2*dim adjacent subdomains
    nProcsPerDim::Vector{Int}   # length of MPI grid in each dimension
    gridPosition::Vector{Int}   # 1-based position of my domain
end

function factorize(n::Int)
    n <= 0 && throw(ArgumentError("must be a positive integer"))
    n == 1 && return [1]
    j = 1
    primeFactors = zeros(Int, Int(ceil(log2(abs(n))))+1)

    for i in 2:Int(ceil(sqrt(n)))
        while n%i==0
            primeFactors[j] = i
            n = div(n, i)
            j+=1
        end
    end
    if n!=1
        primeFactors[j] = n
    else
        j -= 1
    end
    return primeFactors[1:j]
end

function bestCombination(factors::Vector, boxLengths::Vector)
    nDim = size(boxLengths)
    nDim == 1 && return [prod(factors)]

    newFactors = ones(Int, size(factors)+nDim-1)
    newFactors[1:size(factors)] .= factors

    bestCombo = ones(Int, nDim)
    bestCombo[nDim] = prod(factors)
    bestSurfaceArea = Inf64

    for partition in partitions(newFactors, nDim)
        combo = prod.(partition)
        localBoxLengths = boxLengths ./ combo
        halfSurfaceArea = sum(prod(localBoxLengths) ./ localBoxLengths)
        if halfSurfaceArea < bestSurfaceArea
            bestCombo = combo
            bestSurfaceArea = halfSurfaceArea
        end
    end
    return bestCombo
end

function setupSubdomain(me::Int, nProcs::Int, comm::MPI.Comm, boxBounds::Matrix{Float64})
    boxLengths = boxBounds[:, 2] - boxBounds[:, 1]
    nDim = size(boxLengths)
    primeFactors = factorize(nProcs)
    nProcsPerDim = bestCombination(primeFactors, boxLengths)
    procGridPos = ones(Int, nDim)
    allGridPos = ones(Int, nDim, nProcs)

    for i=1:subdomain.nProcs
        allGridPos[:, i] = procGridPos

        procGridPos[nDim] += 1
        for d=nDim:-1:2
            procGridPos[d] <= nProcsPerDim[d] && break

            procGridPos[d] = 1
            procGridPos[d-1] += 1
        end
    end

    gridPosition = allGridPos[:, me+1]
    neighborProcs = zeros(Int, nDim, 2)

    for (i, ind) in enumerate(eachcol(allGridPos))
        for d=1:nDim
            indLo = indHi = ind
            indLo[d] = indLo[d]==1 ? nProcsPerDim[d] : indLo[d]-1
            indHi[d] = indHi[d]==nProcsPerDim[d] ? 1 : indLo[d]+1
            if gridPosition == indLo
                neighborProcs[d, 2] = i
            end
            if gridPosition == indHi
                neighborProcs[d, 1] = i
            end
        end
    end

    localBounds = zeros(Rational, nDim, 2)
    for d in 1:nDim
        localBounds[d, 1] = (gridPosition[d]-1) // nProcsPerDim[d]
        localBounds[d, 2] = gridPosition[d] // nProcsPerDim[d]
    end

    return Subdomain(
        me,
        nProcs,
        comm,
        localBounds,
        boxBounds,
        neighborProcs,
        nProcsPerDim,
        gridPosition
    )
end
