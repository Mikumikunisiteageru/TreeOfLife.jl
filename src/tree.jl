# src/tree.jl

export gettips, gettipnames, alldistinct
export getage, getages
export getparent, getchildren
export preorder, postorder
export isroot, istip, getname, hassibling
export rename, rename!
export getsubtree, getmrca, ismonophyl, getphylodiv
export cutfromroot, cutfromtips
export isbinary
export isisomorph
export getdescs, getdescnames

# CALIBRATING BRANCH LENGTHS (FOR CHRONOTREES)

"""
	calibrate_t_root!(tree::ChronoTree) :: ChronoTree
	calibrate_t_root!(tree::AbstractTree) :: AbstractTree

Calculate all `t_root` values according to `t_branch` values. Used in 
[`fromnewick`](@ref).
"""
function calibrate_t_root!(tree::ChronoTree)
	for node = tree[2:end]
		node.t_root = node.t_branch + tree[node.i_parent].t_root
	end
	tree
end
calibrate_t_root!(tree::AbstractTree) = tree

"""
	calibrate_t_branch!(tree::ChronoTree) :: ChronoTree
	calibrate_t_branch!(tree::AbstractTree) :: AbstractTree

Calibrate all `t_root` values of nodes of the tree so that the root's `t_root` 
is zero, and then recalculate all `t_branch` values according to the new 
`t_root` values. Used in [`getsubtree`](@ref).
"""
function calibrate_t_branch!(tree::ChronoTree)
	t_root_offset = tree[1].t_root
	tree[1].t_root = 0.0
	tree[1].t_branch = 0.0
	for i = preorder(tree)[2:end]
		tree[i].t_root -= t_root_offset
		tree[i].t_branch = tree[i].t_root - 
			tree[tree[i].i_parent].t_root
	end
	tree
end
calibrate_t_branch!(tree::AbstractTree) = tree

# NODE TESTS

"""
	isroot(tree::AbstractTree, i::Integer) :: Bool

Test if the `i`-th node of the tree is the root.
"""
isroot(tree::AbstractTree, i::Integer) = tree[i].i_parent == 0

"""
	istip(tree::AbstractTree, i::Integer) :: Bool

Test if the `i`-th node of the tree is a tip or leaf node.
"""
istip(tree::AbstractTree, i::Integer) = tree[i].i_child == 0

"""
	hassibling(tree::AbstractTree, i::Integer) :: Bool

Test if the `i`-th node of the tree has following sibling(s).
"""
hassibling(tree::AbstractTree, i::Integer) = tree[i].i_sibling > 0

"""
	getname(node::AbstractNode) :: String
	getname(tree::AbstractTree, i::Integer) :: String

Extract the name of the given node as a string.
"""
getname(node::AbstractNode) = node.name
getname(tree::AbstractTree, i::Integer) = getname(tree[i])

"""
	getparent(tree::AbstractTree, i::Integer) :: Int

Find the parent or direct ancestor of `tree[i]` and return its index.
"""
getparent(tree::AbstractTree, i::Integer) = tree[i].i_parent

"""
	getchildren(tree::AbstractTree, i::Integer) :: Vector{Int}

Find all children or direct descendents of `tree[i]` and return their indices.
"""
function getchildren(tree::AbstractTree, i::Integer)
	children = Int[]
	child = tree[i].i_child
	while true
		child == 0 ? (return children) : push!(children, child)
		child = tree[child].i_sibling
	end
end

# GETTING TIPS AND AGE

"""
	gettips(tree::AbstractTree) :: Vector{Int}

Return the indices of all tip nodes of the tree.
"""
gettips(tree::AbstractTree) = filter(i -> istip(tree, i), eachindex(tree))

"""
	gettipnames(tree::AbstractTree) :: Vector{String}

Return the names of all tip nodes of the tree.
"""
gettipnames(tree::AbstractTree) = getname.(tree[gettips(tree)])

"""
	alldistinct(tree::AbstractTree) :: Bool

Test if all tip nodes of the tree have distinct names.
"""
alldistinct(tree::AbstractTree) = allunique(gettipnames(tree))

"""
	isbinary(tree::AbstractTree) :: Bool

Test if the tree is strictly binary or dichotonous, i.e., all non-tip nodes 
have exactly two descendents. 
"""
isbinary(tree::AbstractTree) = all(get_counts(tree, gettipnames(tree)) .== 2)

"""
	mean_(a::Vector{Float64}) :: Float64

Compute the arithmetic mean of a vector of 64-bit float numbers.
"""
mean_(a::Vector{Float64}) = sum(a) / length(a)

"""
	getage(tree::ChronoTree; 
		average=mean_, getrelerr::Bool=false, reltol=1e-8) :: Float64

Return an average age (from the root node) of tip nodes of the tree. 

The argument `average` defines the method for summarizing the ages to one; by 
default it is set to `mean_`.

The argument `getrelerr` controls whether the relative standard deviation is 
appended to the output (`(mean, relstd)`) or not (only `mean`); by default it 
is set to `false`.

The argument `reltol` is a tolerance of relative error. By default it is set 
to `1e-8`. To suppress the judgment, set `reltol=NaN`.
"""
function getage(tree::ChronoTree; 
		average=mean_, getrelerr::Bool=false, reltol=1e-8)
	ages = getfield.(tree[gettips(tree)], :t_root)
	mean = average(ages)
	relerr = maximum(abs.(ages .- mean)) / mean
	! isnan(reltol) && relerr > reltol && throw(ArgumentError(
		"Ages of the tips have relative error $relerr, considered different!"))
	if getrelerr
		return mean, relerr
	else
		# @info "Relative error is $relstd."
		return mean
	end
end

"""
	getages(tree::ChronoTree; average=mean_, reltol=1e-8) :: Vector{Float64}

Return ages (from the root node) of all nodes of the tree.

The argument `average` defines the method for summarizing the ages to one; by 
default it is set to `mean_`.

The argument `reltol` is a tolerance of relative error. By default it is set 
to `1e-8`. To suppress the judgment, set `reltol=NaN`.
"""
function getages(tree::ChronoTree; average=mean_, reltol=1e-8)
	ages = getfield.(tree[gettips(tree)], :t_root)
	mean = average(ages)
	relerr = maximum(abs.(ages .- mean)) / mean
	! isnan(reltol) && relerr > reltol && throw(ArgumentError(
		"Ages of the tips have relative error $relerr, considered different!"))
	return mean .- getfield.(tree, :t_root)
end

# TREE TRAVERSALS

"""
	preorder(tree::AbstractTree, i=1) :: Vector{Int}

Return the pre-order traversal sequence of the whole tree, or its subtree 
with root node `tree[i]`.
"""
function preorder(tree::AbstractTree, i=1)
	sequence = [i]
	tree[i].i_child > 0 && preorder!(sequence, tree, tree[i].i_child)
	sequence
end

"""
	preorder!(sequence, tree::AbstractTree, i=1) :: Nothing

Append the pre-order traversal sequence of the whole tree, or its subtree 
with root node `tree[i]`.
"""
function preorder!(sequence, tree::AbstractTree, i=1) :: Nothing
	push!(sequence, i)
	tree[i].i_child > 0 && preorder!(sequence, tree, tree[i].i_child)
	tree[i].i_sibling > 0 && preorder!(sequence, tree, tree[i].i_sibling)
	nothing
end

"""
	postorder(tree::AbstractTree, i=1) :: Vector{Int}

Return the post-order traversal sequence of the whole tree, or its subtree 
with root node `tree[i]`.
"""
function postorder(tree::AbstractTree, i=1)
	sequence = Int[]
	tree[i].i_child > 0 && postorder!(sequence, tree, tree[i].i_child)
	push!(sequence, i)
end

"""
	postorder!(sequence, tree::AbstractTree, i=1) :: Nothing

Append the post-order traversal sequence of the whole tree, or its subtree 
with root node `tree[i]`.
"""
function postorder!(sequence, tree::AbstractTree, i=1) :: Nothing
	tree[i].i_child > 0 && postorder!(sequence, tree, tree[i].i_child)
	push!(sequence, i)
	tree[i].i_sibling > 0 && postorder!(sequence, tree, tree[i].i_sibling)
	nothing
end

# RENAMING TREE NODES (FOR NEXUS FORMAT)

"""
	rename(oldtree::AbstractTree, 
		oldtonew::Dict{<:AbstractString,<:AbstractString}

Create a new tree whose nodes are respectively renamed from the `oldtree` by 
a mapping from old names to new names. Specifically, nodes with empty names 
remain.
"""
function rename(oldtree::AbstractTree, 
		oldtonew::Dict{<:AbstractString,<:AbstractString})
	newtree = deepcopy(oldtree)
	for node = newtree
		isempty(node.name) && continue
		node.name = oldtonew[node.name]
	end
	newtree
end

"""
	rename!(tree::Tree, 
		oldtonew::Dict{<:AbstractString,<:AbstractString}

Rename nodes of the tree in place by a mapping from old names to new names. 
Specifically, nodes with empty names remain.
"""
function rename!(tree::AbstractTree, 
		oldtonew::Dict{<:AbstractString,<:AbstractString})
	for node = tree
		isempty(node.name) && continue
		node.name = oldtonew[node.name]
	end
	tree
end

# GETTING TEMPORAL SECTIONS (FOR CHRONOGRAMS)

"""
	cutfromroot(tree::ChronoTree, dist::Real; keep::Symbol=:both)
		:: Union{Vector{NTuple{2,Int}}, Vector{Int}}

Find the temporal section by `dist` time units after the root is born. 
The argument `keep` can be set among three options, i.e., `:both` (tuples 
containing parents and childs), `:parent`, and `:child`. 
"""
function cutfromroot(tree::ChronoTree, dist::Real; keep::Symbol=:both)
	keep in [:both, :parent, :child] ||
		throw(ArgumentError(
			"The argument `keep` has to be `:both`, `:parent`, or `:child`!"))
	pcpairs = NTuple{2,Int}[]
	for i = 2:length(tree)
		node = tree[i]
		tree[node.i_parent].t_root < dist <= node.t_root && 
			push!(pcpairs, (node.i_parent, i))
	end
	keep == :parent && return unique(first.(pcpairs))
	keep == :child  && return last.(pcpairs)
	return pcpairs
end

"""
	cutfromtips(tree::ChronoTree, dist::Real; keep::Symbol=:both)
		:: Union{Vector{NTuple{2,Int}}, Vector{Int}}

Find the temporal section by `dist` time units before the root is born.
The argument `keep` can be set among three options, i.e., `:both` (tuples 
containing parents and childs), `:parent`, and `:child`. 
"""
function cutfromtips(tree::ChronoTree, dist::Real; keep::Symbol=:both)
	age = getage(tree; average=minimum)
	cutfromroot(tree, age - dist; keep=keep)
end

# ABOUT TIP SUBSETS

"""
	get_selected(oldtree::AbstractTree, tipset; 
		simplify::Bool=true, keeproot::Bool=false) :: Vector{Int}

Select nodes of a subtree generated from a given set of tips of the tree. Used 
in [`getsubtree`](@ref) and [`isbinary`](@ref). 

Arguments `simplify` and `keeproot` have same meanings as in [`getsubtree`](@ref). 
"""
function get_counts(oldtree::AbstractTree, tipset)
	counts = zeros(Int, length(oldtree))
	for i = postorder(oldtree)[1:end-1]
		if oldtree[i].i_child == 0 && oldtree[i].name in tipset
			counts[i] = 2
		end
		counts[i] >= 1 && (counts[oldtree[i].i_parent] += 1)
	end
	return counts
end

"""
	getsubtree(oldtree::AbstractTree, tipset; 
		simplify::Bool=true, keeproot::Bool=false) :: AbstractTree

Extract the subtree generated from a given set of tips of the tree. 

The argument `simplify` controls whether internal node with only one child 
needs to be reduced, i.e., connecting directly its child and its parent; by 
default it is set to `true`. 

The argument `keeproot` controls whether the original root node needs to be 
contained in the subtree; by default it is set to `false`, in other words, 
yielding a truly minimum spanning tree (MST). 

When `simplify` is set to `false`, the value of `keeproot` has no effect.
"""
function getsubtree(oldtree::AbstractTree, tipset; 
		simplify::Bool=true, keeproot::Bool=false)
	counts = get_counts(oldtree, tipset)
	keeproot && (counts[1] = 2)
	selected = counts .>= simplify + 1
	newnum = sum(selected)
	oldtonew = Dict(findall(selected) .=> 1:newnum)
	newtree = empty(oldtree)
	lastchild = zeros(Int, newnum)
	p = 0
	for oldnode = oldtree[selected]
		p += 1
		newnode = deepcopy(oldnode)
		push!(newtree, newnode)
		newnode.i_sibling = 0
		newnode.i_child = 0
		i = oldnode.i_parent
		while i > 0 && ! selected[i]
			i = oldtree[i].i_parent
		end
		if i == 0
			newnode.i_parent = 0
			continue
		end
		newnode.i_parent = oldtonew[i]
		if newtree[newnode.i_parent].i_child == 0
			newtree[newnode.i_parent].i_child = p
		end
		if lastchild[newnode.i_parent] > 0
			newtree[lastchild[newnode.i_parent]].i_sibling = p
		end
		lastchild[newnode.i_parent] = p
	end
	return calibrate_t_branch!(newtree)
end

"""
	getmrca(tree::AbstractTree, tipset) :: Int

Find the index of the most recent common ancestor node for a set of nodes. 
"""
getmrca(tree::AbstractTree, tipset) = findfirst(get_counts(tree, tipset) .>= 2)

"""
	getdescs(tree::AbstractTree, mrca::Int) :: Vector{Int}

Find the indices of all tip descendents from a common ancestor node.
"""
getdescs(tree::AbstractTree, mrca::Int) = 
	filter(i -> istip(tree, i), preorder(tree, mrca))

"""
	getdescnames(tree::AbstractTree, mrca::Int) :: Vector{String}
	getdescnames(tree::AbstractTree) :: Vector{Vector{String}}

Find the names of all tip descendents from a common ancestor node, or such 
descendent name lists for all nodes of the tree.
"""
getdescnames(tree::AbstractTree, mrca::Int) = 
	getname.(tree[getdescs(tree, mrca)])
function getdescnames(tree::AbstractTree)
	descnames = [String[] for _ = tree]
	for i = postorder(tree)[1:end-1]
		istip(tree, i) && push!(descnames[i], getname(tree, i))
		append!(descnames[tree[i].i_parent], descnames[i])
	end
	descnames
end

"""
	ismonophyl(tree::AbstractTree, tipset) :: Bool

Test if a given set of tip nodes are monophyletic based on the tree.
"""
ismonophyl(tree::AbstractTree, tipset) = 
	Set(getdescnames(tree, getmrca(tree, tipset))) == Set(tipset)

"""
	sum_t_branch(tree::ChronoTree)

Compute the sum of branch lengths of the tree. Used in [`getphylodiv`](@ref). 
"""
sum_t_branch(tree::ChronoTree) = sum(n.t_branch for n = tree)

"""
	getphylodiv(tree::ChronoTree, tipset; keeproot::Bool=false)

Compute the phylogenetic diversity (PD) of a given set of tips of the tree, 
i.e., the sum of branch lengths of the subtree generated from the set. 

The argument `keeproot` controls whether the original root node needs to be 
contained in the subtree; by default it is set to `false`.
"""
getphylodiv(tree::ChronoTree, tipset; keeproot::Bool=false) = 
	sum_t_branch(getsubtree(tree, tipset, simplify=true, keeproot=keeproot))

# ISOMORPHISM

"""
	tree_hash(tree::CladoTree, h::UInt=zero(UInt)) :: UInt
	tree_hash(tree::ChronoTree, h::UInt=zero(UInt)) :: UInt

Compute a hash value for a phylogenetic tree so that isomorphic trees 
necessarily have the same hash value (tested by [`isisomorph`](@ref)).
"""
function tree_hash(tree::CladoTree, h::UInt=zero(UInt))
	hashes = fill(UInt(1), length(tree))
	for i = eachindex(tree)[end:-1:2]
		hashes[i] *= hash(tree[i].name, h)
		hashes[tree[i].i_parent] += hashes[i]
	end
	return hashes[1] * hash(tree[1].name, h)
end
function tree_hash(tree::ChronoTree, h::UInt=zero(UInt))
	hashes = fill(UInt(1), length(tree))
	for i = eachindex(tree)[end:-1:2]
		hashes[i] *= xor(hash(tree[i].name, h), hash(tree[i].t_branch, h)) 
		hashes[tree[i].i_parent] += hashes[i]
	end
	return hashes[1] * hash(tree[1].name, h)
end

"""
	isisomorph(tree1::CladoTree, tree2::CladoTree) :: Bool
	isisomorph(tree1::ChronoTree, tree2::ChronoTree) :: Bool

Test if two trees are isomorphic. 

When both phylogenetic tree are dated, the isomorphism implies that branch 
lengths are correspondingly equal; otherwise, only the tree topology are 
compared.
"""
isisomorph(tree1::CladoTree, tree2::CladoTree) = 
	tree_hash(tree1) == tree_hash(tree2)
isisomorph(tree1::ChronoTree, tree2::ChronoTree) = 
	tree_hash(tree1) == tree_hash(tree2)
isisomorph(tree1::AbstractTree, tree2::AbstractTree) = 
	tree_hash(CladoTree(tree1)) == tree_hash(CladoTree(tree2))
