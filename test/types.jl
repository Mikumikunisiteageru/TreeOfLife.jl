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
	@test (tree[1] = ChronoNode(name="A")) == ChronoNode(name="A")
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
	@test CladoTree(CladoTree()) == CladoTree()
	@test CladoTree([CladoNode(name="A")]) == CladoTree([CladoNode(name="A")])
	@test ChronoTree([ChronoNode(name="A")]) == 
		ChronoTree([ChronoNode(name="A")])
	@test CladoNode(ChronoNode(name="A")) == CladoNode(name="A")
	@test CladoTree(tree) == 
		CladoTree([CladoNode(name="A"), CladoNode(name="B")])
	@test_throws KeyError tree[:flowercolor]
	@test_throws MethodError tree[:color] = Set(:color)
	@test (tree[1, :color] = "blue-green") == "blue-green"
	@test tree[1, :color] == "blue-green"
	@test tree[:color] == Dict(1 => "blue-green")
	@test tree.traitvalues == Dict((1, :color) => "blue-green")
	@test tree.traits == Set([:color])
	@test haskey(tree, :flowercolor) == false
	@test haskey(tree, 1, :flowercolor) == false
	@test haskey(tree, :color) == true
	@test haskey(tree, 1, :color) == true
	@test haskey(tree, 2, :color) == false
end
