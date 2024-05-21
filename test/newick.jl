# test/newick.jl

using TreeOfLife
const TOL = TreeOfLife
using Test

@testset "ignore_comments" begin
	@test TOL.ignore_comments("1[2]3") == "13"
	@test TOL.ignore_comments("[2]3") == "3"
	@test TOL.ignore_comments("1[2]") == "1"
	@test TOL.ignore_comments("[2]") == ""
	@test TOL.ignore_comments("[1[2]3]") == ""
	@test TOL.ignore_comments("生命[之]树") == "生命树"
end

@testset "from_newick" begin
	@test_throws ArgumentError TOL.from_newick("(A,B,(C,D))")
	@test TOL.from_newick("(A,B,(C,D));") ==
		['(', "A", "B", '(', "C", "D", ')', ')']
	@test TOL.from_newick("(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);") ==
		['(', "A", 0.1, "B", 0.2, '(', "C", 0.3, "D", 0.4, ')', 0.5, ')']
	@test TOL.from_newick("(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;") ==
		['(', "A", 0.1, "B", 0.2, 
			  '(', "C", 0.3, "D", 0.4, ')', ("E",), 0.5, ')', ("F",)]
	@test TOL.from_newick("(A:0.1,B:0.1):0.0;") ==
		['(', "A", 0.1, "B", 0.1, ')']
	@test TOL.from_newick("(A:0.1,B:0.1):0.2;") ==
		['(', "A", 0.1, "B", 0.1, ')', 0.2]
end

@testset "fromnewick" begin
	global tree_b = fromnewick("(A,B,(C,D)E)F;")
	@test isa(tree_b, CladoTree)
	@test_throws ArgumentError fromnewick("(A:0.1,B,(C,D));")
	global tree = fromnewick("(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;")
	@test isa(tree, ChronoTree)
	@test length(tree) == 6
	@test getfield.(tree, :name) == ["F", "A", "B", "E", "C", "D"]
	moraceae = String(read("files/moraceae.tre"))
	global tree1 = fromnewick(moraceae)
	@test length(tree1) == 639
	@test tree1[1].name == ""
	@test tree1[2].name == "[&rate=0.0010896447096906676]"
	@test tree1[3].name == "296[&rate=9.19731143726293E-4]"
	global tree2 = fromnewick(moraceae, nocomments=true)
	@test length(tree2) == 639
	@test tree2[1].name == ""
	@test tree2[2].name == ""
	@test tree2[3].name == "296"
	@test isa(fromnewick("(A:0.1,B:0.1):0.0;"), ChronoTree)
	@test_throws ArgumentError fromnewick("(A:0.1,B:0.1):0.2;")
end

@testset "tonewick" begin
	@test tonewick(fromnewick("(A,B,(C,D));")) == "(A,B,(C,D));"
	@test tonewick(tree_b) == "(A,B,(C,D)E)F;"
	@test tonewick(tree) == "(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;"
end

@testset "read/write newick" begin
	filename = joinpath(tempdir(), "test.tre")
	writenewick(filename, tree)
	newtree = readnewick(filename)
	@test tree == newtree
	writenewick(filename, tree_b)
	newtree_b = readnewick(filename)
	@test tree_b == newtree_b
end
