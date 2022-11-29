module TreeOfLife

export Node, Tree, ChronoNode, ChronoTree
export fromnewick, gettips, getage
export preorder, postorder
export rename, rename!
export readnexus
export subtree, phylodiv

using Parameters

abstract type Node end
abstract type Tree end

@with_kw mutable struct ChronoNode <: Node
	name::String = ""
	i_parent::Int = 0
	i_sibling::Int = 0
	i_child::Int = 0
	t_root::Float64 = 0.0
	t_branch::Float64 = 0.0
end

struct ChronoTree <: Tree
	nodes::Vector{ChronoNode}
end
ChronoTree() = ChronoTree(Vector{ChronoNode}())

Base.length(tree::Tree) = length(tree.nodes)
Base.getindex(tree::Tree, i) = getindex(tree.nodes, i)
Base.setindex!(tree::Tree, node::Node, i) = setindex!(tree.nodes, node, i)
Base.firstindex(tree::Tree) = firstindex(tree.nodes)
Base.lastindex(tree::Tree) = lastindex(tree.nodes)
Base.iterate(tree::Tree) = iterate(tree.nodes)
Base.iterate(tree::Tree, i) = iterate(tree.nodes, i)
Base.eachindex(tree::Tree) = eachindex(tree.nodes)
Base.push!(tree::Tree, node::Node) = push!(tree.nodes, node)

"""
	Base.:(==)(node1::ChronoNode, node2::ChronoNode)

Tells if two nodes are identical, in the sense that they have the same name, 
the same parent, sibling, and child, as well as they have approximate branch 
lengths.
"""
function Base.:(==)(node1::ChronoNode, node2::ChronoNode)
	node1.name      == node2.name      &&
	node1.i_parent  == node2.i_parent  &&
	node1.i_sibling == node2.i_sibling &&
	node1.i_child   == node2.i_child   &&
	isapprox(node1.t_root,   node2.t_root) &&
	isapprox(node1.t_branch, node2.t_branch)
end

"""
	Base.:(==)(tree1::ChronoTree, tree2::ChronoTree)

Tells if two trees are identical, in the sense that they have the same nodes. 
Cf. `isisomorph`.
"""
Base.:(==)(tree1::ChronoTree, tree2::ChronoTree) = 
	length(tree1) == length(tree2) && all(tree1 .== tree2)

"""
	ignore_comments(str::AbstractString)

Gets rid of comments enclosed by square brackets (possibly nested) from `str`.
"""
function ignore_comments(str::AbstractString)
	output = ""
	depth = 0
	for char = str
		char == '[' && (depth += 1)
		depth == 0 && (output *= char)
		char == ']' && (depth -= 1)
	end
	@assert depth == 0
	return output
end

"""
	parse_newick(str::AbstractString)

Parses a newick-format tree string into segments of different types according 
to their syntactical meanings. Used by `fromnewick`.
"""
function parse_newick(str::AbstractString)
	elements = Union{Char, String, Float64, Tuple{String}}[]
	register = ""
	for char = strip(str)
		if occursin(char, "(:,);")
			if ! isempty(register)
				push!(elements, register)
				register = ""
			end
			push!(elements, char)
		else
			register *= char
		end
	end
	for i = findall(elements .== ':') .+ 1
		elements[i] = parse(Float64, elements[i])
	end
	for i = findall(elements .== ')') .+ 1
		if isa(elements[i], String)
			elements[i] = (elements[i],)
		end
	end
	@assert elements[end] == ';'
	filter(x -> ! in(x, [':', ',', ';']), elements)
end

"""
	fromnewick(str::AbstractString; nocomments::Bool=false)

Converts a newick-format tree string `str` to a `ChronoTree` instance.
The argument `nocomments` controls whether the comments (enclosed by square 
brackets) are wiped out; by default it is set to `false`, i.e., all comments 
are kept.
"""
function fromnewick(str::AbstractString; nocomments::Bool=false)
	tree = ChronoTree([ChronoNode()])
	stack = [1, 0]
	p = 1
	elements = parse_newick(str)
	decorate = nocomments ? ignore_comments : identity
	for element = elements[2:end]
		if isa(element, String) || element == '('
			p += 1
			if stack[end] == 0
				tree[stack[end-1]].i_child = p
			else
				tree[stack[end]].i_sibling = p
			end
			push!(tree, ChronoNode(i_parent=stack[end-1]))
			stack[end] = p
			if isa(element, String)
				tree[end].name = decorate(element)
			else
				push!(stack, 0)
			end
		elseif element == ')'
			pop!(stack)
		elseif isa(element, Tuple{String})
			tree[stack[end]].name = decorate(element[1])
		else # if isa(element, Float64)
			tree[stack[end]].t_branch = element
		end
	end
	for node = tree[2:end]
		node.t_root = node.t_branch + tree[node.i_parent].t_root
	end
	tree
end

"""
	gettips(tree::ChronoTree)

Returns indices of tips nodes on `tree`.
"""
gettips(tree::ChronoTree) = filter(i -> tree[i].i_child == 0, eachindex(tree))

"""
	getage(tree::ChronoTree)

Returns an arithmetic mean of ages (from the root node) of tip nodes on `tree`.
The argument `getrelstd` controls whether the relative standard deviation is 
appended to the output (`(mean, relstd)`) or not (only `mean`); by default it 
is set to `false`.
"""
function getage(tree::ChronoTree; getrelstd::Bool=false)
	ages = getfield.(gettips(tree), :t_root)
	mean = sum(ages) / length(ages)
	std = sqrt(sum((ages .- mean) .^ 2) / (length(ages) - 1))
	relstd = std / mean
	if getrelstd
		return mean, relstd
	else
		@info "Relative standard deviation is $relstd."
		return mean
	end
end

"""
	preorder(tree::ChronoTree, i=1)

Returns the pre-order traversal sequence of the whole `tree`, or its subtree 
with root node `tree[i]`.
"""
function preorder(tree::ChronoTree, i=1)
	sequence = [i]
	tree[i].i_child > 0 && preorder!(sequence, tree, tree[i].i_child)
	sequence
end
function preorder!(sequence, tree::ChronoTree, i=1)
	push!(sequence, i)
	tree[i].i_child > 0 && preorder!(sequence, tree, tree[i].i_child)
	tree[i].i_sibling > 0 && preorder!(sequence, tree, tree[i].i_sibling)
end

"""
	preorder(tree::ChronoTree, i=1)

Returns the post-order traversal sequence of the whole `tree`, or its subtree 
with root node `tree[i]`.
"""
function postorder(tree::ChronoTree, i=1)
	sequence = Int[]
	tree[i].i_child > 0 && postorder!(sequence, tree, tree[i].i_child)
	push!(sequence, i)
end
function postorder!(sequence, tree::ChronoTree, i=1)
	tree[i].i_child > 0 && postorder!(sequence, tree, tree[i].i_child)
	push!(sequence, i)
	tree[i].i_sibling > 0 && postorder!(sequence, tree, tree[i].i_sibling)
end

"""
	rename(oldtree::ChronoTree, 
		oldtonew::Dict{<:AbstractString,<:AbstractString}

Creates a new tree whose nodes are respectively renamed from the `oldtree` by 
a mapping from old names to new names. Specifically, nodes with empty names 
remain.
"""
function rename(oldtree::ChronoTree, 
		oldtonew::Dict{<:AbstractString,<:AbstractString})
	newtree = deepcopy(oldtree)
	for node = newtree
		isempty(node.name) && continue
		node.name = oldtonew[node.name]
	end
	newtree
end

"""
	rename!(tree::ChronoTree, 
		oldtonew::Dict{<:AbstractString,<:AbstractString}

Renames nodes in `tree` by a mapping from old names to new names. 
Specifically, nodes with empty names remain.
"""
function rename!(tree::ChronoTree, 
		oldtonew::Dict{<:AbstractString,<:AbstractString})
	for node = tree
		isempty(node.name) && continue
		node.name = oldtonew[node.name]
	end
	tree
end

function readnexus(fname::AbstractString; every=0)
	trees = ChronoTree[]
	open(fname, "r") do fin
		line = readline(fin)
		@assert line == "#NEXUS"
		oldtonew = Dict{String,String}()
		stage = 0
		treecnt = 0
		while ! eof(fin)
			line = readline(fin)
			if line == "\tTranslate"
				stage = 1
				continue
			elseif stage == 1 && line == "\t\t;"
				stage = 2
				continue
			elseif stage == 2 && line == "End;"
				stage = 3
				break
			end
			if stage == 1
				old, new = split(strip(line, ['\t', ',']), ' ')
				oldtonew[old] = new
			elseif stage == 2
				str = last(split(strip(line, '\t'), ' '))
				tree = fromnewick(str, nocomments=true)
				push!(trees, rename!(tree, oldtonew))
				every > 0 && (treecnt += 1) % every == 0 && println(treecnt)
			end
		end
	end
	return trees
end

"""
	get_selected(oldtree::ChronoTree, tipset::Set{<:AbstractString};
		simplify::Bool=true, keeproot::Bool=false)

Selects nodes of a subtree generated from a given set of tips on `tree`. Used by `subtree`. 
Arguments `simplify` and `keeproot` have same meanings as in `subtree`. 
"""
function get_selected(oldtree::ChronoTree, tipset::Set{<:AbstractString};
		simplify::Bool=true, keeproot::Bool=false)
	counts = zeros(Int, length(oldtree))
	for i = postorder(oldtree)[1:end-1]
		if oldtree[i].i_child == 0 && oldtree[i].name in tipset
			counts[i] = 2
		end
		counts[i] >= 1 && (counts[oldtree[i].i_parent] += 1)
	end
	keeproot && (counts[1] = 2)
	counts .>= simplify + 1
end

"""
	calibrate_time!(newtree::ChronoTree)

Calibrates all `t_root` values of nodes in `newtree`, so that the root's 
`t_root` is zero. Then, all `t_branch` values are recalculated according to 
`t_root` values. Used by `subtree`.
"""
function calibrate_time!(newtree::ChronoTree)
	t_root_offset = newtree[1].t_root
	newtree[1].t_root = 0.0
	newtree[1].t_branch = 0.0
	for i = preorder(newtree)[2:end]
		newtree[i].t_root -= t_root_offset
		newtree[i].t_branch = newtree[i].t_root - 
			newtree[newtree[i].i_parent].t_root
	end
	newtree
end

"""
	subtree(oldtree::ChronoTree, tipset::Set{<:AbstractString}; 
		simplify::Bool=true, keeproot::Bool=false)

Extracts the subtree generated from a given set of tips on `tree`. 
The argument `simplify` controls whether internal node with only one child needs to be reduced, i.e., connecting directly its child and its parent; by default it is set to `true`. 
The argument `keeproot` controls whether the original root node needs to be contained in the subtree; by default it is set to `false`, in other words, yielding a truly minimum spanning tree (MST). 
When `simplify` is set to `false`, the value of `keeproot` has no effect.
"""
function subtree(oldtree::ChronoTree, tipset::Set{<:AbstractString}; 
		simplify::Bool=true, keeproot::Bool=false)
	selected = get_selected(oldtree, tipset; 
		simplify=simplify, keeproot=keeproot)
	newnum = sum(selected)
	oldtonew = Dict(findall(selected) .=> 1:newnum)
	newtree = ChronoTree()
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
	calibrate_time!(newtree)
end

"""
	sum_t_branch(tree::ChronoTree)

Calculates the sum of branch lengths of `tree`. Used by `phylodiv`. 
"""
sum_t_branch(tree::ChronoTree) = sum(n.t_branch for n = tree)

"""
	phylodiv(tree::ChronoTree, tipset::Set{<:AbstractString}; 
		keeproot::Bool=false)

Calculates the phylogenetic diversity (PD) of a given set of tips on `tree`, 
i.e., the sum of branch lengths of the subtree generated from the set. 
The argument `keeproot` controls whether the original root node needs to be 
contained in the subtree; by default it is set to `false`.
"""
phylodiv(tree::ChronoTree, tipset::Set{<:AbstractString}; 
		keeproot::Bool=false) = 
	sum_t_branch(subtree(tree, tipset, simplify=true, keeproot=keeproot))

end # module TreeOfLife
