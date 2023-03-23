# src/types.jl

using Parameters

export AbstractNode, ChronoNode, CladoNode
export AbstractTree, ChronoTree, CladoTree

abstract type AbstractNode end

abstract type AbstractTree end

@with_kw mutable struct ChronoNode <: AbstractNode
	name::String = ""
	i_parent::Int = 0
	i_sibling::Int = 0
	i_child::Int = 0
	t_root::Float64 = 0.0
	t_branch::Float64 = 0.0
end

@with_kw mutable struct CladoNode <: AbstractNode
	name::String = ""
	i_parent::Int = 0
	i_sibling::Int = 0
	i_child::Int = 0
end
CladoNode(node::ChronoNode) = CladoNode(
	name=node.name, i_parent=node.i_parent, 
	i_sibling=node.i_sibling, i_child=node.i_child)

struct ChronoTree <: AbstractTree
	nodes::Vector{ChronoNode}
end
ChronoTree() = ChronoTree(Vector{ChronoNode}())

struct CladoTree <: AbstractTree
	nodes::Vector{CladoNode}
end
CladoTree() = CladoTree(Vector{CladoNode}())
CladoTree(tree::ChronoTree) = CladoTree(CladoNode.(tree))

Base.length(tree::AbstractTree) = length(tree.nodes)
Base.getindex(tree::AbstractTree, i) = getindex(tree.nodes, i)
Base.setindex!(tree::AbstractTree, node::AbstractNode, i) = 
	setindex!(tree.nodes, node, i)
Base.firstindex(tree::AbstractTree) = firstindex(tree.nodes)
Base.lastindex(tree::AbstractTree) = lastindex(tree.nodes)
Base.iterate(tree::AbstractTree) = iterate(tree.nodes)
Base.iterate(tree::AbstractTree, i) = iterate(tree.nodes, i)
Base.eachindex(tree::AbstractTree) = eachindex(tree.nodes)
Base.push!(tree::AbstractTree, node::AbstractNode) = push!(tree.nodes, node)
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
	==(tree1::ChronoTree, tree2::ChronoTree) :: Bool

Test if two trees are identical, in the sense that they have the same nodes. 

Identical trees are always isomorphic (tested by [`isomorphic`](@ref)).
"""
Base.:(==)(tree1::AbstractTree, tree2::AbstractTree) = 
	length(tree1) == length(tree2) && all(tree1 .== tree2)
