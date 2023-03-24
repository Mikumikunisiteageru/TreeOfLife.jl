# src/types.jl

using Parameters

export AbstractNode, ChronoNode, CladoNode
export AbstractTree, ChronoTree, CladoTree

"""
	AbstractNode

Abstract supertype for [`ChronoNode`](@ref) and [`CladoNode`](@ref).
"""
abstract type AbstractNode end

"""
	AbstractTree

Abstract supertype for [`ChronoTree`](@ref) and [`CladoTree`](@ref).
"""
abstract type AbstractTree end

"""
	ChronoNode <: AbstractNode

Type for nodes in a [`ChronoTree`](@ref). Compare [`CladoNode`](@ref).
"""
@with_kw mutable struct ChronoNode <: AbstractNode
	name::String = ""
	i_parent::Int = 0
	i_sibling::Int = 0
	i_child::Int = 0
	t_root::Float64 = 0.0
	t_branch::Float64 = 0.0
end

"""
	CladoNode <: AbstractNode
	
	CladoNode(node::ChronoNode) :: CladoNode

Type for nodes in a [`CladoTree`](@ref). 

A [`ChronoNode`](@ref) can be converted to a `CladoNode` by removing all 
information about time.
"""
@with_kw mutable struct CladoNode <: AbstractNode
	name::String = ""
	i_parent::Int = 0
	i_sibling::Int = 0
	i_child::Int = 0
end
CladoNode(node::ChronoNode) = CladoNode(
	name=node.name, i_parent=node.i_parent, 
	i_sibling=node.i_sibling, i_child=node.i_child)

"""
	ChronoTree <: AbstractTree
	
	ChronoTree(nodes::AbstractVector{ChronoNode}) :: ChronoTree
	ChronoTree() :: ChronoTree

Type for chronograms or dated phylogenetic trees, assumed to be rooted, 
comprising [`ChronoNode`](@ref) instances. Compare [`CladoTree`](@ref). 

When called with no arguments, the constructor returns an empty `ChronoTree`.
"""
struct ChronoTree <: AbstractTree
	nodes::Vector{ChronoNode}
	traits::Set{Symbol}
	traitvalues::Dict{Tuple{Int, Symbol}, Any}
	ChronoTree(nodes::AbstractVector{ChronoNode}) = 
		new(Vector(nodes), Set{Symbol}(), Dict{Tuple{Int, Symbol}, Any}())
end
ChronoTree() = ChronoTree(Vector{ChronoNode}())

"""
	CladoTree <: AbstractTree
	
	CladoTree() :: CladoTree
	CladoTree(nodes::AbstractVector{CladoNode}) :: CladoTree
	CladoTree(tree::CladoTree) :: CladoTree
	CladoTree(tree::ChronoTree) :: CladoTree

Type for cladograms or undated phylogenetic trees, assumed to be rooted, 
comprising [`CladoNode`](@ref) instances. 

When called with no arguments, the constructor returns an empty `CladoTree`.

A [`ChronoTree`](@ref) can be converted to a `CladoTree` by removing all 
information about time.
"""
struct CladoTree <: AbstractTree
	nodes::Vector{CladoNode}
	traits::Set{Symbol}
	traitvalues::Dict{Tuple{Int, Symbol}, Any}
	CladoTree(nodes::AbstractVector{CladoNode}) = 
		new(Vector(nodes), Set{Symbol}(), Dict{Tuple{Int, Symbol}, Any}())
end
CladoTree() = CladoTree(Vector{CladoNode}())
CladoTree(tree::CladoTree) = tree
CladoTree(tree::ChronoTree) = CladoTree(CladoNode.(tree))

"""
	length(tree::AbstractTree) :: Int

Return the number of nodes in a phylogenetic tree.
"""
Base.length(tree::AbstractTree) = length(tree.nodes)

"""
	getindex(tree::AbstractTree, i::Integer) :: AbstractNode
	getindex(tree::AbstractTree, i::Integer, trait::Symbol) :: Any
	getindex(tree::AbstractTree, trait::Symbol) :: Dict{Int, <:Any}

Get the `i`-th node of the tree, or its trait, or a `Dict` containing the 
trait of all nodes of the tree. See also [`setindex!`](@ref).
"""
Base.getindex(tree::AbstractTree, i::Integer) = getindex(tree.nodes, i)
Base.getindex(tree::AbstractTree, i::Integer, trait::Symbol) = 
	getindex(tree.traitvalues, i, trait)
Base.getindex(tree::AbstractTree, trait::Symbol) = 
	trait in tree.traits ? 
		Dict(k[1] => v for (k, v) = tree.traitvalues if k[2] == trait) : 
		throw(KeyError(trait))

"""
	setindex!(tree::AbstractTree, node::AbstractNode, i::Integer) :: Vector
	setindex!(tree::AbstractTree, tv, i::Integer, trait::Symbol) :: Dict

Set the `i`-th node of the tree or its trait. See also [`getindex`](@ref).
"""
Base.setindex!(tree::AbstractTree, node::AbstractNode, i::Integer) = 
	setindex!(tree.nodes, node, i)
function Base.setindex!(tree::AbstractTree, tv, i::Integer, trait::Symbol)
	push!(tree.traits, trait)
	setindex!(tree.traitvalues, tv, i, trait)
end

Base.firstindex(tree::AbstractTree) = firstindex(tree.nodes)
Base.lastindex(tree::AbstractTree) = lastindex(tree.nodes)
Base.iterate(tree::AbstractTree) = iterate(tree.nodes)
Base.iterate(tree::AbstractTree, i) = iterate(tree.nodes, i)
Base.eachindex(tree::AbstractTree) = eachindex(tree.nodes)
Base.push!(tree::AbstractTree, node::AbstractNode) = push!(tree.nodes, node)

"""
	empty(tree::AbstractTree) :: AbstractTree

Construct a phylogenetic tree with only the root node of the same type.
"""
Base.empty(tree::AbstractTree) = typeof(tree)()

"""
	==(node1::CladoNode, node2::CladoNode) :: Bool
	==(node1::ChronoNode, node2::ChronoNode) :: Bool

Test if two nodes are identical, in the sense that they have the same name, 
the same parent, sibling, and child, as well as they have approximate branch 
lengths.
"""
function Base.:(==)(node1::CladoNode, node2::CladoNode)
	node1.name      == node2.name      &&
	node1.i_parent  == node2.i_parent  &&
	node1.i_sibling == node2.i_sibling &&
	node1.i_child   == node2.i_child
end
function Base.:(==)(node1::ChronoNode, node2::ChronoNode)
	node1.name      == node2.name      &&
	node1.i_parent  == node2.i_parent  &&
	node1.i_sibling == node2.i_sibling &&
	node1.i_child   == node2.i_child   &&
	isapprox(node1.t_root,   node2.t_root) &&
	isapprox(node1.t_branch, node2.t_branch)
end

"""
	==(tree1::AbstractTree, tree2::AbstractTree) :: Bool

Test if two trees are identical, in the sense that they are of the same type, 
isomorphic, and have the same node ordering; specifically, for dated trees or 
[`ChronoTree`](@ref) instances, the node times are correspondingly equal. 

Identical trees are always isomorphic (can be tested by [`isomorphic`](@ref)).
"""
Base.:(==)(tree1::AbstractTree, tree2::AbstractTree) = 
	typeof(tree1) == typeof(tree2) && 
	length(tree1) == length(tree2) && all(tree1 .== tree2)
