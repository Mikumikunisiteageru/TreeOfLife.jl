# src/tree.jl

export fromnewick, tonewick
export gettips, gettipnames, alldistinct
export getage, getages
export preorder, postorder
export isroot, istip, getname, hassibling
export rename, rename!
export subtree, getmrca, ismonophyl, phylodiv
export cutfromroot, cutfromtips
export isbinary
export treehash, isomorphic
export getdescnames

# CALIBRATING BRANCH LENGTHS (FOR CHRONOTREES)

"""
	calibrate_t_root!(tree::ChronoTree)
	calibrate_t_root!(tree::Tree)

Calculates all `t_root` values according to `t_branch` values. Used by 
`fromnewick`.
"""
function calibrate_t_root!(tree::ChronoTree)
	for node = tree[2:end]
		node.t_root = node.t_branch + tree[node.i_parent].t_root
	end
	tree
end

calibrate_t_root!(tree::AbstractTree) = tree

"""
	calibrate_t_branch!(tree::ChronoTree)
	calibrate_t_branch!(tree::Tree)

Calibrates all `t_root` values of nodes in `tree`, so that the root's `t_root` 
is zero, and then recalculates all `t_branch` values according to the new 
`t_root` values. Used by `subtree`.
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

# PARSING NEWICK-FORMAT TREE STRINGS

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
	from_newick(str::AbstractString)

Parses a Newick-format tree string into segments of different types according 
to their syntactical meanings. Used by `fromnewick`.
"""
function from_newick(str::AbstractString)
	endswith(str, ';') || 
		throw(ArgumentError(
			"The Newick-format string does not end with a semicolon!"))
	elements = [] # Union{Char, String, Float64, Tuple{String}}
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
	get_tree_type(elements::Vector)

Tells whether a vector of elements parsed from some Newick-syntax tree string 
is a chronogram or a cladogram. Used by `fromnewick`.
"""
function get_tree_type(elements::Vector)
	n_left  = sum(elements .== '(')
	n_right = sum(elements .== ')')
	n_left == n_right || 
		throw(ArgumentError(
			"Numbers of left and right parentheses are unequal!"))
	n_tip    = sum(isa.(elements, String))
	n_branch = sum(isa.(elements, Float64))
	n_branch == 0 && 
		return CladoTree,  CladoNode
	n_branch == n_tip + n_left - 1 && 
		return ChronoTree, ChronoNode
	throw(ArgumentError(
		"Branch lengths should occur on either all or none of the nodes!"))
end

"""
	fromnewick(str::AbstractString; nocomments::Bool=false)

Converts a Newick-format tree string `str` to a `ChronoTree` instance.
The argument `nocomments` controls whether the comments (enclosed by square 
brackets) are wiped out; by default it is set to `false`, i.e., all comments 
are kept.
"""
function fromnewick(str::AbstractString; nocomments::Bool=false)
	elements = from_newick(str)
	T, N = get_tree_type(elements)
	tree = T([N()])
	stack = [1, 0]
	p = 1
	decorate = nocomments ? ignore_comments : identity
	for element = elements[2:end]
		if isa(element, String) || element == '('
			p += 1
			if stack[end] == 0
				tree[stack[end-1]].i_child = p
			else
				tree[stack[end]].i_sibling = p
			end
			push!(tree, N(i_parent=stack[end-1]))
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
	return calibrate_t_root!(tree)
end

function to_newick!(elements::Vector, tree::CladoTree, i::Int)
	ic = tree[i].i_child
	if ic == 0
		push!(elements, tree[i].name)
	else # ic > 0
		push!(elements, '(')
		to_newick!(elements, tree, ic)
		push!(elements, ')', tree[i].name)
	end
	is = tree[i].i_sibling
	if is > 0
		push!(elements, ',')
		to_newick!(elements, tree, is)
	end
end

function to_newick!(elements::Vector, tree::ChronoTree, i::Int)
	ic = tree[i].i_child
	if ic == 0
		push!(elements, tree[i].name, ':', tree[i].t_branch)
	else # ic > 0
		push!(elements, '(')
		to_newick!(elements, tree, ic)
		push!(elements, ')', tree[i].name, ':', tree[i].t_branch)
	end
	is = tree[i].i_sibling
	if is > 0
		push!(elements, ',')
		to_newick!(elements, tree, is)
	end
end

function tonewick(tree::CladoTree)
	elements = []
	to_newick!(elements, tree, 1)
	push!(elements, ';')
	return join(elements)
end

function tonewick(tree::ChronoTree)
	elements = []
	to_newick!(elements, tree, 1)
	pop!(elements)
	elements[end] = ';'
	return join(elements)
end

# NODE JUDGMENTS

isroot(tree::AbstractTree, i::Int) = i == 1
istip(tree::AbstractTree, i::Int) = tree[i].i_child == 0
hassibling(tree::AbstractTree, i::Int) = tree[i].i_sibling > 0
getname(tree::AbstractTree, i::Int) = tree[i].name
getname(node::AbstractNode) = node.name

# GETTING TIPS AND AGE

"""
	gettips(tree::Tree)

Returns indices of tips nodes on `tree`.
"""
gettips(tree::AbstractTree) = 
	filter(i -> tree[i].i_child == 0, eachindex(tree))

gettipnames(tree::AbstractTree) = getname.(tree[gettips(tree)])

alldistinct(tree::AbstractTree) = allunique(gettipnames(tree))

isbinary(tree::AbstractTree) = all(get_counts(tree, gettipnames(tree)) .== 2)

"""
	mean(a::Vector{Float64})

Calculates the arithmetic mean of a vector of float numbers.
"""
mean(a::Vector{Float64}) = sum(a) / length(a)

"""
	getage(tree::ChronoTree; 
		average=mean, getrelerr::Bool=false, reltol=1e-8)

Returns an average age (from the root node) of tip nodes on `tree`. 
The argument `average` defines the method for summarizing the ages to one; by 
default it is set to `mean`.
The argument `getrelerr` controls whether the relative standard deviation is 
appended to the output (`(mean, relstd)`) or not (only `mean`); by default it 
is set to `false`.
The argument `reltol` is a tolerance of relative error. By default it is set 
to `1e-8`. To suppress the judgment, set `reltol=NaN`.
"""
function getage(tree::ChronoTree; 
		average=mean, getrelerr::Bool=false, reltol=1e-8)
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

function getages(tree::ChronoTree; average=mean, reltol=1e-8)
	ages = getfield.(tree[gettips(tree)], :t_root)
	mean = average(ages)
	relerr = maximum(abs.(ages .- mean)) / mean
	! isnan(reltol) && relerr > reltol && throw(ArgumentError(
		"Ages of the tips have relative error $relerr, considered different!"))
	return mean .- getfield.(tree, :t_root)
end

# TREE TRAVERSALS

"""
	preorder(tree::Tree, i=1)

Returns the pre-order traversal sequence of the whole `tree`, or its subtree 
with root node `tree[i]`.
"""
function preorder(tree::AbstractTree, i=1)
	sequence = [i]
	tree[i].i_child > 0 && preorder!(sequence, tree, tree[i].i_child)
	sequence
end
function preorder!(sequence, tree::AbstractTree, i=1)
	push!(sequence, i)
	tree[i].i_child > 0 && preorder!(sequence, tree, tree[i].i_child)
	tree[i].i_sibling > 0 && preorder!(sequence, tree, tree[i].i_sibling)
end

"""
	preorder(tree::Tree, i=1)

Returns the post-order traversal sequence of the whole `tree`, or its subtree 
with root node `tree[i]`.
"""
function postorder(tree::AbstractTree, i=1)
	sequence = Int[]
	tree[i].i_child > 0 && postorder!(sequence, tree, tree[i].i_child)
	push!(sequence, i)
end
function postorder!(sequence, tree::AbstractTree, i=1)
	tree[i].i_child > 0 && postorder!(sequence, tree, tree[i].i_child)
	push!(sequence, i)
	tree[i].i_sibling > 0 && postorder!(sequence, tree, tree[i].i_sibling)
end

# RENAMING TREE NODES (FOR NEXUS FORMAT)

"""
	rename(oldtree::Tree, 
		oldtonew::Dict{<:AbstractString,<:AbstractString}

Creates a new tree whose nodes are respectively renamed from the `oldtree` by 
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

Renames nodes in `tree` by a mapping from old names to new names. 
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

Finds the temporal section by `dist` units after the root is born. 
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

Finds the temporal section by `dist` units before the first tip is born.
The argument `keep` can be set among three options, i.e., `:both` (tuples 
containing parents and childs), `:parent`, and `:child`. 
"""
function cutfromtips(tree::ChronoTree, dist::Real; keep::Symbol=:both)
	age = getage(tree; average=minimum)
	cutfromroot(tree, age - dist; keep=keep)
end

# ABOUT TIP SUBSETS

"""
	get_selected(oldtree::Tree, tipset; 
		simplify::Bool=true, keeproot::Bool=false)

Selects nodes of a subtree generated from a given set of tips on `tree`. Used 
by `subtree`. 
Arguments `simplify` and `keeproot` have same meanings as in `subtree`. 
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
	subtree(oldtree::ChronoTree, tipset; 
		simplify::Bool=true, keeproot::Bool=false)

Extracts the subtree generated from a given set of tips on `tree`. 
The argument `simplify` controls whether internal node with only one child 
needs to be reduced, i.e., connecting directly its child and its parent; by 
default it is set to `true`. 
The argument `keeproot` controls whether the original root node needs to be 
contained in the subtree; by default it is set to `false`, in other words, 
yielding a truly minimum spanning tree (MST). 
When `simplify` is set to `false`, the value of `keeproot` has no effect.
"""
function subtree(oldtree::AbstractTree, tipset; 
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

getmrca(tree::AbstractTree, tipset) = findfirst(get_counts(tree, tipset) .>= 2)

getdescs(tree::AbstractTree, mrca::Int) = 
	filter(i -> istip(tree, i), preorder(tree, mrca))

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

ismonophyl(tree::AbstractTree, tipset) = 
	Set(getdescnames(tree, getmrca(tree, tipset))) == Set(tipset)

"""
	sum_t_branch(tree::ChronoTree)

Calculates the sum of branch lengths of `tree`. Used by `phylodiv`. 
"""
sum_t_branch(tree::ChronoTree) = sum(n.t_branch for n = tree)

"""
	phylodiv(tree::ChronoTree, tipset; keeproot::Bool=false)

Calculates the phylogenetic diversity (PD) of a given set of tips on `tree`, 
i.e., the sum of branch lengths of the subtree generated from the set. 
The argument `keeproot` controls whether the original root node needs to be 
contained in the subtree; by default it is set to `false`.
"""
phylodiv(tree::ChronoTree, tipset; keeproot::Bool=false) = 
	sum_t_branch(subtree(tree, tipset, simplify=true, keeproot=keeproot))

# ISOMORPHISM

function treehash(tree::CladoTree, h::UInt=zero(UInt))
	hashes = fill(UInt(1), length(tree))
	for i = eachindex(tree)[end:-1:2]
		hashes[i] *= hash(tree[i].name, h)
		hashes[tree[i].i_parent] += hashes[i]
	end
	return hashes[1] * hash(tree[1].name, h)
end

function treehash(tree::ChronoTree, h::UInt=zero(UInt))
	hashes = fill(UInt(1), length(tree))
	for i = eachindex(tree)[end:-1:2]
		hashes[i] *= xor(hash(tree[i].name, h), hash(tree[i].t_branch, h)) 
		hashes[tree[i].i_parent] += hashes[i]
	end
	return hashes[1] * hash(tree[1].name, h)
end

isomorphic(tree1::CladoTree, tree2::CladoTree) = 
	treehash(tree1) == treehash(tree2)

isomorphic(tree1::ChronoTree, tree2::ChronoTree) = 
	treehash(tree1) == treehash(tree2)
