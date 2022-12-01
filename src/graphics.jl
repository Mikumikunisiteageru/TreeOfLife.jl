# src/graphics.jl

export setplot, showtree

struct PlotModule end

function plot_unset(a...; kw...)
	throw(ErrorException(
		"Plot module unspecified! Use e.g. `using PyPlot; setplot(PyPlot)`."))
end

Base.getproperty(::PlotModule, ::Symbol) = plot_unset

setplot(x::Module) = global Plot = x

Plot = PlotModule()
plot(a...; kw...) = Plot.plot(a...; kw...)
text(a...; kw...) = Plot.text(a...; kw...)
gca(a...; kw...) = Plot.gca(a...; kw...)
axis(a...; kw...) = Plot.axis(a...; kw...)

function calc_pos(tree::CladoTree; x0=0.0, y0=0.0, levelsep=1, siblingsep=-1)
	levelsep <= 0 && throw(DomainError(levelsep, 
		"The argument `levelsep` has to be positive!"))
	siblingsep == 0 && throw(DomainError(siblingsep, 
		"The argument `siblingsep` has to be non-zero!"))
	n = length(tree)
	x = fill(x0, n)
	y = Vector{Float64}(undef, n)
	yc0 = fill(+Inf, n)
	yc1 = fill(-Inf, n)
	y1 = y0
	for i = postorder(tree)
		if istip(tree, i)
			y[i] = y1 += siblingsep
		else
			x[i] -= levelsep
			y[i] = (yc0[i] + yc1[i]) / 2
		end
		(p = tree[i].i_parent) == 0 && continue
		x[p] = min(x[p], x[i])
		yc0[p] = min(yc0[p], y[i])
		yc1[p] = max(yc1[p], y[i])
	end
	return x, y
end

function claim_length(array::AbstractArray, name::String, n::Int)
	length(array) == n || throw(ArgumentError(
		"The length of `$name` is different from that of `tree`!"))
	return i -> array[i]
end

claim_length(element, name, n) = i -> element

function draw_tree(tree::CladoTree, x::Vector{Float64}, y::Vector{Float64};
		linewidth=0.8, linecolor="k")
	n = length(tree)
	getlw = claim_length(linewidth, "linewidth", n)
	getlc = claim_length(linecolor, "linecolor", n)
	for i = preorder(tree)[2:end]
		p = tree[i].i_parent
		plot(x[[p,p,i]], y[[p,i,i]], "-"; color=getlc(i), lw=getlw(i))
	end
end

function draw_text(tree::CladoTree, x::Vector{Float64}, y::Vector{Float64}; 
		textsep=0.0, fontsize="small", fontcolor="k")
	for i = eachindex(tree)
		istip(tree, i) || continue
		text(x[i]+textsep, y[i], "  " * tree[i].name; 
			fontsize=fontsize, color=fontcolor, ha="left", va="center")
	end
end

function showtree(tree::CladoTree; x0=0.0, y0=0.0, levelsep=1, siblingsep=-1, 
		linewidth=0.8, linecolor="k", 
		textsep=0.0, fontsize="small", fontcolor="k")
	x, y = calc_pos(tree; x0=x0, y0=y0, 
		levelsep=levelsep, siblingsep=siblingsep)
	draw_tree(tree, x, y; linewidth=linewidth, linecolor=linecolor)
	draw_text(tree, x, y; 
		textsep=textsep, fontsize=fontsize, fontcolor=fontcolor)
	gca().set_position([0, 0, 1, 1])
	axis("off")
end
