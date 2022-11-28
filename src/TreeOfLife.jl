module TreeOfLife

export ChronoNode, Chronogram
export fromnewick, gettips, getage
export preorder, postorder
export rename, rename!
export readnexus, filter

using Parameters

@with_kw mutable struct ChronoNode
	name::String = ""
	i_parent::Int = 0
	i_sibling::Int = 0
	i_child::Int = 0
	t_root::Float64 = 0.0
	t_branch::Float64 = 0.0
end
const Chronogram = Vector{ChronoNode}

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

function fromnewick(str::AbstractString; nocomments::Bool=false)
	tree = [ChronoNode()]
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

gettips(tree::Chronogram) = filter(i -> tree[i].i_child == 0, eachindex(tree))

function getage(tree::Chronogram)
	ages = getfield.(gettips(tree), :t_root)
	mean = sum(ages) / length(ages)
	std = sqrt(sum((ages .- mean) .^ 2) / (length(ages) - 1))
	relstd = std / mean
	@info "Relative standard deviation is $relstd."
	return mean
end

"""
	preorder(tree::Chronogram, i=1)

Returns the pre-order traversal sequence of the whole `tree`, or its subtree 
with root node `tree[i]`.
"""
function preorder(tree::Chronogram, i=1)
	sequence = [i]
	tree[i].i_child > 0 && preorder!(sequence, tree, tree[i].i_child)
	sequence
end
function preorder!(sequence, tree::Chronogram, i=1)
	push!(sequence, i)
	tree[i].i_child > 0 && preorder!(sequence, tree, tree[i].i_child)
	tree[i].i_sibling > 0 && preorder!(sequence, tree, tree[i].i_sibling)
end

"""
	preorder(tree::Chronogram, i=1)

Returns the post-order traversal sequence of the whole `tree`, or its subtree 
with root node `tree[i]`.
"""
function postorder(tree::Chronogram, i=1)
	sequence = Int[]
	tree[i].i_child > 0 && postorder!(sequence, tree, tree[i].i_child)
	push!(sequence, i)
end
function postorder!(sequence, tree::Chronogram, i=1)
	tree[i].i_child > 0 && postorder!(sequence, tree, tree[i].i_child)
	push!(sequence, i)
	tree[i].i_sibling > 0 && postorder!(sequence, tree, tree[i].i_sibling)
end

"""
	rename(oldtree::Chronogram, 
		oldtonew::Dict{<:AbstractString,<:AbstractString}

Creates a new tree whose nodes are respectively renamed from the `oldtree` by 
a mapping from old names to new names. Specifically, nodes with empty names 
remain.
"""
function rename(oldtree::Chronogram, 
		oldtonew::Dict{<:AbstractString,<:AbstractString})
	newtree = deepcopy(oldtree)
	for node = newtree
		isempty(node.name) && continue
		node.name = oldtonew[node.name]
	end
	newtree
end

"""
	rename!(tree::Chronogram, 
		oldtonew::Dict{<:AbstractString,<:AbstractString}

Renames nodes in `tree` by a mapping from old names to new names. 
Specifically, nodes with empty names remain.
"""
function rename!(tree::Chronogram, 
		oldtonew::Dict{<:AbstractString,<:AbstractString})
	for node = tree
		isempty(node.name) && continue
		node.name = oldtonew[node.name]
	end
	tree
end

function readnexus(fname::AbstractString; every=0)
	trees = Chronogram[]
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

function filter(tipset::Set{<:AbstractString}, oldtree::Chronogram)
	selected = falses(length(oldtree))
	for i = postorder(oldtree)[1:end-1]
		if oldtree[i].i_child == 0 && oldtree[i].name in tipset
			selected[i] = true
		end
		selected[oldtree[i].i_parent] |= selected[i]
	end
	oldtonew = Dict(findall(selected) .=> 1:sum(selected))
	newtree = Chronogram()
	lastchild = zeros(Int, sum(selected))
	p = 0
	for oldnode = oldtree[selected]
		p += 1
		newnode = deepcopy(oldnode)
		push!(newtree, newnode)
		newnode.i_sibling = 0
		newnode.i_child = 0
		oldnode.i_parent == 0 && continue
		newnode.i_parent = oldtonew[oldnode.i_parent]
		if newtree[newnode.i_parent].i_child == 0
			newtree[newnode.i_parent].i_child = p
		end
		if lastchild[newnode.i_parent] > 0
			newtree[lastchild[newnode.i_parent]].i_sibling = p
		end
		lastchild[newnode.i_parent] = p
	end
	newtree
end

end # module TreeOfLife
