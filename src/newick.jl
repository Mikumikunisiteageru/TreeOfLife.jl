# src/newick.jl

export readnewick, writenewick
export fromnewick, tonewick

"""
	ignore_comments(str::AbstractString) :: String

Remove comments enclosed by square brackets (possibly nested) in a string.
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
	from_newick(str::AbstractString) :: Vector

Parse a Newick-format tree string into segments of different types according 
to their syntactical meanings. Used in [`fromnewick`](@ref).
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
	elements = filter(x -> ! in(x, [':', ',', ';']), elements)
	if isa(elements[end], Float64) && iszero(elements[end])
		return elements[1:end-1]
	else
		return elements
	end
end

"""
	get_tree_type(elements::Vector) :: 
		Tuple{Type{<:AbstractTree}, Type{<:AbstractNode}}

Test whether a vector of elements parsed from some Newick-syntax tree string 
is a chronogram or a cladogram. Used in [`fromnewick`](@ref).
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
	fromnewick(str::AbstractString; nocomments::Bool=false) :: AbstractTree

Create a phylogenetic tree from a Newick-format tree string.

Its inverse function is [`tonewick`](@ref).

The argument `nocomments` controls whether the comments (enclosed by square 
brackets) are wiped out; by default it is set to `false`, i.e., all comments 
are kept.
"""
function fromnewick(str::AbstractString; nocomments::Bool=false)
	elements = from_newick(strip(str))
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

"""
	to_newick!(elements::AbstractVector, 
		tree::CladoTree, i::Integer) :: AbstractVector
	to_newick!(elements::AbstractVector, 
		tree::ChronoTree, i::Integer) :: AbstractVector

Convert a phylogenetic tree to segments of a Newick-format string. Used in 
[`tonewick`](@ref).
"""
function to_newick!(elements::AbstractVector, tree::CladoTree, i::Integer)
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
	return elements
end
function to_newick!(elements::AbstractVector, tree::ChronoTree, i::Integer)
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
	return elements
end

"""
	tonewick(tree::CladoTree) :: String
	tonewick(tree::ChronoTree) :: String

Express a phylogenetic tree as a Newick-format string.

Its inverse function is [`fromnewick`](@ref).
"""
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

"""
	readnewick(filename::AbstractString) :: AbstractTree

Read a phylogenetic tree from disk. Compare [`writenewick`](@ref).
"""
readnewick(filename::AbstractString) = 
	fromnewick(strip(read(filename, String)))

"""
	readnewick(filename::AbstractString, tree::AbstractTree) :: AbstractTree

Write a phylogenetic tree to disk. Compare [`readnewick`](@ref).
"""
writenewick(filename::AbstractString, tree::AbstractTree) = 
	write(filename, tonewick(tree))
