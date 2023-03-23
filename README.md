# TreeOfLife.jl

[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://Mikumikunisiteageru.github.io/TreeOfLife.jl/stable)
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://Mikumikunisiteageru.github.io/TreeOfLife.jl/dev)
[![CI](https://github.com/Mikumikunisiteageru/TreeOfLife.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/Mikumikunisiteageru/TreeOfLife.jl/actions/workflows/CI.yml)
[![Codecov](https://codecov.io/gh/Mikumikunisiteageru/TreeOfLife.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/Mikumikunisiteageru/TreeOfLife.jl)
[![Aqua.jl Quality Assurance](https://img.shields.io/badge/Aquajl-%F0%9F%8C%A2-aqua.svg)](https://github.com/JuliaTesting/Aqua.jl)

TreeOfLife.jl defines data types for cladograms and chronograms in phylogenetics and methods analyzing these trees.

In development. 

Alternative packages include:
- [Phylo.jl](https://github.com/EcoJulia/Phylo.jl) ("in beta", and its outdated prototype [Phylogenetics.jl](https://github.com/SabrinaJaye/Phylogenetics.jl))
- [Phylogenies.jl](https://github.com/BioJulia/Phylogenies.jl) ("in development")

## Examples

The followings are some simple examples illustrating the usage of this package. Please see [the documentation](https://Mikumikunisiteageru.github.io/TreeOfLife.jl/) for more details of TreeOfLife.jl.

```julia
julia> using TreeOfLife

julia> tree = fromnewick("(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;")
ChronoTree(ChronoNode[ChronoNode("F", 0, 0, 2, 0.0, 0.0), ChronoNode("A", 1, 3, 0, 0.1, 0.1), ChronoNode("B", 1, 4, 0, 0.2, 0.2), ChronoNode("E", 1, 0, 5, 0.5, 0.5), ChronoNode("C", 4, 6, 0, 0.8, 0.3), ChronoNode("D", 4, 0, 0, 0.9, 0.4)])

julia> tonewick(tree)
"(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;"

julia> tipnames = gettipnames(tree)
4-element Vector{String}:
 "A"
 "B"
 "C"
 "D"

julia> getmrca(tree, tipnames)
1

julia> phylodiv(tree, tipnames)
1.5

julia> ismonophyl(tree, ["A", "B"])
false

julia> ismonophyl(tree, ["C", "D"])
true
```
