export readDatafile

function removeComment(line::AbstractString)
    strip(split(line, "#")[1])
end

# Currently missing the masses
function readDatafile(datafile::String, nDim::Int)
    open(datafile, "r") do f
        readline(f)  # skip first line
        line = readline(f)  # next line blank

        nAtoms = 0
        imageFlag = false
        boxBounds = zeros(Float64, nDim, 2)

        # Parse header, read until non-blank line without header keyword
        while !eof(f)
            line = readline(f)
            line = removeComment(line)
            (line == "") && continue
            
            if occursin("atoms", line)
                nAtoms = parse(Int, split(line)[1])
            elseif occursin("xlo xhi", line)
                values = split(line)
                boxBounds[1, : ] = [parse(Float64, values[1]), parse(Float64, values[2])]
                for i=2:nDim
                    values = split(removeComment(readline(f)))
                    boxBounds[i, : ] = [parse(Float64, values[1]), parse(Float64, values[2])]
                end
            else
                break
            end
        end

        atoms = Atoms(
            nAtoms, 
            0, 
            zeros(Float64, nDim, nAtoms), 
            zeros(Float64, nDim, nAtoms), 
            zeros(Float64, nDim, nAtoms), 
            zeros(Int, 100, nAtoms),
            zeros(Int, nAtoms),
            ones(Int, nAtoms),
            ones(Float64, nAtoms),
            zeros(Int, nDim, nAtoms),
        )

        while !eof(f)
            if removeComment(line) == ""
                line = readline(f)
                continue
            end
            if occursin("Atoms", line)
                readline(f)  # Skip next line

                line = readline(f)
                splitLine = split(removeComment(line))
                if length(splitLine) == 2 + 2*nDim
                    imageFlag = true
                elseif length(splitLine) == 2 + nDim
                    imageFlag = false
                else
                    throw(DimensionMismatch(
                        "invalid number of columns in Atoms section of data file."
                    ))
                end

                for i = 1:nAtoms
                    atoms.ids[i] = parse(Int, splitLine[1])
                    atoms.types[i] = parse(Int, splitLine[2])
                    atoms.positions[:, i] = parse.(Float64, splitLine[3:2+nDim])
                    imageFlag && atoms.images[:, i] = parse.(Int, splitLine[3+nDim:2+2*nDim])
                    splitline = split(removeComment(readline(f)))
                end
            elseif occursin("Velocities", line)
                readline(f)
                line = split(removeComment(readline(f)))
                for i in 1:nAtoms
                    atoms.velocities[:, i] = parse.(Float64, line)
                    line = split(removeComment(readline(f)))
                end
            end
            line = readline(f)
        end
    end
    return (atoms, boxBounds)
end
