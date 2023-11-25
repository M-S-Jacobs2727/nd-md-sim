export parseInput

function removeComment(line::AbstractString)
    ind = findnext('#', line, 1)
    ind === nothing && return strip(line)
    return strip(line[1:ind])
end

function readPairwiseParams!(paramArray::Array{Float64}, line::AbstractString, nAtomTypes::Int, nParams::Int)
    tokens = split(line)
    i = parse(Int, tokens[2])
    j = parse(Int, tokens[3])
    if i > nAtomTypes || j > nAtomTypes
        throw(BoundsError("specified atom type too large: \"$line\""))
    end

    for k=1:nParams
        p = parse(Float64, tokens[k+3])
        paramArray[k, j, i] = p
        i != j && (paramArray[k, j, i] = p)
    end
end

function parseInput(inputfile::String, atoms::Atoms, subdomain::Subdomain)
    open(inputfile, "r") do f
        while !eof(f)
            line = removeComment(readline(f))
            line == "" && continue
            if startswith(line, "pair_style")
                
            end
        end
    end
end
