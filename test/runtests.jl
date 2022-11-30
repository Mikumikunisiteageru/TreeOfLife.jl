using TreeOfLife
const TOL = TreeOfLife
using Test

cd(joinpath(pkgdir(TreeOfLife), "test", "files"))

@testset "TYPES" begin
	@test CladoNode  <: Node
	@test ChronoNode <: Node
	@test CladoTree  <: Tree
	@test ChronoTree <: Tree
	@test length(CladoTree())  == 0
	@test length(ChronoTree()) == 0
	tree = ChronoTree([ChronoNode(name="A"), ChronoNode(name="B")])
	@test length(tree) == 2
	@test tree[1].name == "A"
	@test firstindex(tree) == 1
	@test lastindex(tree) == 2
	@test eachindex(tree) == Base.OneTo(2)
	@test empty(tree) == ChronoTree()
	@test CladoNode(name="A") == CladoNode(name="A")
	@test ChronoNode(name="A") == ChronoNode(name="A")
	@test CladoTree() == CladoTree()
	@test ChronoTree() == ChronoTree()
	@test CladoTree([CladoNode(name="A")]) == CladoTree([CladoNode(name="A")])
	@test ChronoTree([ChronoNode(name="A")]) == 
		ChronoTree([ChronoNode(name="A")])
	@test CladoNode(ChronoNode(name="A")) == CladoNode(name="A")
	@test CladoTree(tree) == 
		CladoTree([CladoNode(name="A"), CladoNode(name="B")])
end

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
end

@testset "fromnewick" begin
	global tree_b = fromnewick("(A,B,(C,D)E)F;")
	@test isa(tree_b, CladoTree)
	@test_throws ArgumentError fromnewick("(A:0.1,B,(C,D));")
	global tree = fromnewick("(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;")
	@test isa(tree, ChronoTree)
	@test length(tree) == 6
	@test getfield.(tree, :name) == ["F", "A", "B", "E", "C", "D"]
	moraceae = String(read("moraceae.tre"))
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
end

@testset "tonewick" begin
	@test tonewick(fromnewick("(A,B,(C,D));")) == "(A,B,(C,D));"
	@test tonewick(tree_b) == "(A,B,(C,D)E)F;"
	@test tonewick(tree) == "(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;"
end

@testset "gettips" begin
	tree = fromnewick("(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;")
	@test gettips(tree) == [2, 3, 5, 6]
end

@testset "subtree" begin
	@test subtree(tree, Set(["A", "B"])) == fromnewick("(A:0.1,B:0.2)F;")
	@test subtree(tree, Set(["C", "D"])) == fromnewick("(C:0.3,D:0.4)E;")
	@test subtree(tree, Set(["C", "D"]); simplify=false) == 
		fromnewick("((C:0.3,D:0.4)E:0.5)F;")
	@test subtree(tree, Set(["C", "D"]); simplify=true, keeproot=true) == 
		fromnewick("((C:0.3,D:0.4)E:0.5)F;")
	@test subtree(tree, Set(["C"]); simplify=true, keeproot=true) == 
		fromnewick("(C:0.8)F;")
	@test subtree(tree, Set(["A", "C"])) == fromnewick("(A:0.1,C:0.8)F;")
	@test subtree(tree, Set(["A", "C"]); simplify=false) == 
		fromnewick("(A:0.1,(C:0.3)E:0.5)F;")
end

@testset "getmrca" begin
	@test getmrca(tree, ["A", "B"]) == 1
	@test getmrca(tree, ["C", "D"]) == 4
	@test getmrca(tree, ["C"]) == 5
	@test getmrca(tree, ["A", "C"]) == 1
	@test isnothing(getmrca(tree, []))
	@test isnothing(getmrca(tree, ["G"]))
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