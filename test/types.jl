# test/types.jl

using TreeOfLife
const TOL = TreeOfLife
using Test

@testset "types" begin
	@test ChronoNode <: AbstractNode
	@test CladoNode  <: AbstractNode
	@test ChronoTree <: AbstractTree
	@test CladoTree  <: AbstractTree
	@test length(ChronoTree()) == 0
	@test length(CladoTree())  == 0
	tree = ChronoTree([ChronoNode(name="A"), ChronoNode(name="B")])
	@test length(tree) == 2
	@test tree[1].name == "A"
	@test firstindex(tree) == 1
	@test lastindex(tree) == 2
	@test eachindex(tree) == Base.OneTo(2)
	@test empty(tree) == ChronoTree()
	@test CladoNode(name="A") == CladoNode(name="A")
	@test ChronoNode(name="A") == ChronoNode(name="A")
	@test CladoTree() == CladoTree()
	@test ChronoTree() == ChronoTree()
	@test CladoNode() != ChronoNode()
	@test CladoTree() != ChronoTree()
	@test CladoTree([CladoNode(name="A")]) == CladoTree([CladoNode(name="A")])
	@test ChronoTree([ChronoNode(name="A")]) == 
		ChronoTree([ChronoNode(name="A")])
	@test CladoNode(ChronoNode(name="A")) == CladoNode(name="A")
	@test CladoTree(tree) == 
		CladoTree([CladoNode(name="A"), CladoNode(name="B")])
end
