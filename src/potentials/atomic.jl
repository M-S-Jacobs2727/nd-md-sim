export AtomicPotential

abstract type AtomicPotential end
# include("atomic/ljcut.jl")
for (root, dirs, files) in walkdir("atomic")
    if files.size == 0
        for file in files
            include(joinpath(root, file))
        end
    end
end
