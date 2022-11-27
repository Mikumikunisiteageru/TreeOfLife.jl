module TreeOfLife

export ChronoNode, Chronogram
export read_newick

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

function parse_newick(str::AbstractString)
	elements = Union{Char, String, Float64}[]
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
	@assert elements[[1, end-1, end]] == ['(', ')', ';']
	filter(x -> ! in(x, [':', ',']), elements[2:end-2])
end

function update_time!(tree::Chronogram)
	for node = tree[2:end]
		node.t_root = node.t_branch + tree[node.i_parent].t_root
	end
	tree
end

function read_newick(str::AbstractString)
	tree = [ChronoNode()]
	stack = [[1], Int[]]
	p = 1
	elements = parse_newick(str)
	for element = elements
		if isa(element, String) || element == '('
			p += 1
			if isempty(stack[end])
				tree[stack[end-1][end]].i_child = p
			else
				tree[stack[end][end]].i_sibling = p
			end
			push!(tree, ChronoNode(i_parent=stack[end-1][end]))
			push!(stack[end], p)
			if isa(element, String)
				tree[end].name = element
			else
				push!(stack, Int[])
			end
		elseif element == ')'
			pop!(stack)
		elseif isa(element, Float64)
			tree[stack[end][end]].t_branch = element
		end
	end
	update_time!(tree)
end

end # module TreeOfLife
