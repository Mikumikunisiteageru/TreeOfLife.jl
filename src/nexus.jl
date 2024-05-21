# src/nexus.jl

export readnexus

"""
	readnexus(filename::AbstractString; every=0) :: Vector{ChronoTree}

Read dated phylogenetic trees from a Nexus-format file on disk. 

!!! warning
	Under development; correctness not guaranteed.
"""
function readnexus(filename::AbstractString; every=0)
	trees = ChronoTree[]
	open(filename, "r") do fin
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
			elseif stage == 1 && strip(line) == ";"
				stage = 2
				continue
			elseif stage == 2 && line == "End;"
				stage = 3
				break
			end
			if stage == 1
				old, new = split(strip(line, ['\t', ',', ' ']), ' ')
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
