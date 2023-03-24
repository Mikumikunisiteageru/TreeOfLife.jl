# test/tree.jl

using TreeOfLife
const TOL = TreeOfLife
using Test

tree = fromnewick("(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;")

@testset "gettips" begin
	@test gettips(tree) == [2, 3, 5, 6]
end

@testset "getparent" begin
	@test getparent.([tree], 1:6) == [0, 1, 1, 1, 4, 4]
end

@testset "getchildren" begin
	@test getchildren.([tree], 1:6) == [[2,3,4], [], [], [5,6], [], []]
end

@testset "getsubtree" begin
	@test getsubtree(tree, Set(["A", "B"])) == fromnewick("(A:0.1,B:0.2)F;")
	@test getsubtree(tree, Set(["C", "D"])) == fromnewick("(C:0.3,D:0.4)E;")
	@test getsubtree(tree, Set(["C", "D"]); simplify=false) == 
		fromnewick("((C:0.3,D:0.4)E:0.5)F;")
	@test getsubtree(tree, Set(["C", "D"]); simplify=true, keeproot=true) == 
		fromnewick("((C:0.3,D:0.4)E:0.5)F;")
	@test getsubtree(tree, Set(["C"]); simplify=true, keeproot=true) == 
		fromnewick("(C:0.8)F;")
	@test getsubtree(tree, Set(["A", "C"])) == fromnewick("(A:0.1,C:0.8)F;")
	@test getsubtree(tree, Set(["A", "C"]); simplify=false) == 
		fromnewick("(A:0.1,(C:0.3)E:0.5)F;")
end

@testset "getmrca" begin
	@test getmrca(tree, ["A", "B"]) == 1
	@test getmrca(tree, ["C", "D"]) == 4
	@test getmrca(tree, ["C"]) == 5
	@test getmrca(tree, ["A", "C"]) == 1
	@test isa(getmrca(tree, []), Nothing)
	@test isa(getmrca(tree, ["G"]), Nothing)
end

@testset "getdescnames" begin
	@test getdescnames(tree) == 
		[["A", "B", "C", "D"], ["A"], ["B"], ["C", "D"], ["C"], ["D"]]
	@test getdescnames(tree) == getdescnames.([tree], eachindex(tree))
end

@testset "ismonophyl" begin
	@test ismonophyl(tree, ["A", "B"]) == false
	@test ismonophyl(tree, ["C", "D"]) == true
	@test ismonophyl(tree, ["A", "C"]) == false
	@test ismonophyl(tree, ["C"]) == true
end

@testset "cutfromroot" begin
	@test cutfromroot(tree, 0.0) == []
	@test cutfromroot(tree, 0.05) == [(1,2), (1,3), (1,4)]
	@test cutfromroot(tree, 0.1) == [(1,2), (1,3), (1,4)]
	@test cutfromroot(tree, 0.15) == [(1,3), (1,4)]
	@test cutfromroot(tree, 0.05, keep=:both) == [(1,2), (1,3), (1,4)]
	@test cutfromroot(tree, 0.05, keep=:parent) == [1]
	@test cutfromroot(tree, 0.05, keep=:child) == [2, 3, 4]
end