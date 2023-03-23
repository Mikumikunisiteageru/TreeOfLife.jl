# src/TreeOfLife.jl

module TreeOfLife
@doc let
    path = joinpath(dirname(@__DIR__), "README.md")
    include_dependency(path)
    read(path, String)
end TreeOfLife

include("types.jl")

include("newick.jl")
include("nexus.jl")

include("tree.jl")
include("trees.jl")

# include("graphics.jl")

end # module TreeOfLife
