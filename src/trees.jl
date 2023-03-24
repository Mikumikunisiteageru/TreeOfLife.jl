# src/trees.jl

export getconsensus

"""
	count_clade_ages(trees::Vector{ChronoTree}; every::Int=0) 
		:: Dict{Set{String}, Vector{Float64}}

Count the age for every possible tip node combination as long as occured at 
least once in some of the trees. 

The argument `every`, if set to a positive integer, makes the (probably long) 
analyzing process print a log per `every` trees; the default value is `0`.
"""
function count_clade_ages(trees::Vector{ChronoTree}; every::Int=0)
	counter = Dict{Set{String},Vector{Float64}}()
	ntree = length(trees)
	for i = 1:ntree
		every > 0 && i % every == 0 && println(i, '/', ntree)
		for (d, a) = zip(getdescnames(trees[i]), getages(trees[i]))
			s = Set(d)
			push!(get!(counter, s, Float64[]), a)
		end
	end
	counter
end

"""
	count_clades(trees::Vector{<:AbstractTree}; every::Int=0)
		:: Dict{Set{String}, Int}

Count every possible tip node combination as long as occured at least once in 
some of the trees. Used in [`getconsensus`](@ref).

The argument `every`, if set to a positive integer, makes the (probably long) 
analyzing process print a log per `every` trees; the default value is `0`.
"""
function count_clades(trees::Vector{<:AbstractTree}; every::Int=0)
	counter = Dict{Set{String},Int}()
	ntree = length(trees)
	for i = 1:ntree
		every > 0 && i % every == 0 && println(i, '/', ntree)
		for d = getdescnames(trees[i])
			s = Set(d)
			counter[s] = get(counter, s, 0) + 1
		end
	end
	counter
end

"""
	construct_tree(parents::Vector{Int}) :: Tuple{CladoTree, Int}

Construct a cladogram from parent indices. Used in [`getconsensus`](@ref).
"""
function construct_tree(parents::Vector{Int})
	m = length(parents)
	roots = findall(parents .== 0)
	length(roots) == 1 || throw(ArgumentError(
		"The graph is disconnected and thus not a rooted tree!"))
	root = roots[1]
	tree = CladoTree([CladoNode(i_parent=parents[i]) for i = 1:m])
	for i = 1:m
		(p = parents[i]) == 0 && continue
		tree[i].i_sibling = tree[p].i_child
		tree[p].i_child = i
	end
	return tree, root
end

"""
	getconsensus(trees::Vector{<:AbstractTree}; 
		threshold::Float64=0.5, every::Int=0) :: CladoTree

Summarize phylogenetic trees to one such that clades with support rate no less 
than a threshold remain in the consensus tree. 

The argument `threshold` should be a number between `0.5` and `1.0`; the default 
value is `0.5`.

The argument `every`, if set to a positive integer, makes the (probably long) 
analyzing process print a log per `every` trees; the default value is `0`.
"""
function getconsensus(trees::Vector{<:AbstractTree}; 
		threshold::Float64=0.5, every::Int=0)
	0.5 <= threshold <= 1.0 || throw(ArgumentError(
		"The argument `threshold` has to be in [0.5, 1.0]!"))
	0.5 == threshold && (threshold += eps(0.5))
	nthres = ceil(Int, threshold * length(trees))
	clades = Set{String}[]
	counter = count_clades(trees; every=every)
	for s = keys(counter)
		counter[s] >= nthres && push!(clades, s)
	end
	sort!(clades, by=length)
	clades_len = length.(clades)
	mclades = length(clades)
	parents = fill(0, mclades)
	for i = 1 : mclades-1
		for j = i+1 : mclades
			clades_len[i] == clades_len[j] && continue
			if issubset(clades[i], clades[j])
				parents[i] = j
				break
			end
		end
	end
	oldtree, root = construct_tree(parents)
	for i = 1:mclades
		length(clades[i]) == 1 && (oldtree[i].name = first(clades[i]))
	end
	newtoold = preorder(oldtree, root)
	oldtonew = Dict(0 => 0)
	setindex!.([oldtonew], 1:mclades, newtoold);
	CladoTree([CladoNode(name=n.name, 
		i_parent=oldtonew[n.i_parent], 
		i_sibling=oldtonew[n.i_sibling], 
		i_child=oldtonew[n.i_child]) for n = oldtree[newtoold]])
end
