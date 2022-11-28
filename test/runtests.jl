using TreeOfLife
const TOL = TreeOfLife
using Test

cd(joinpath(pkgdir(TreeOfLife), "test", "files"))

@testset "ignore_comments" begin
	@test TOL.ignore_comments("1[2]3") == "13"
	@test TOL.ignore_comments("[2]3") == "3"
	@test TOL.ignore_comments("1[2]") == "1"
	@test TOL.ignore_comments("[2]") == ""
	@test TOL.ignore_comments("[1[2]3]") == ""
	@test TOL.ignore_comments("生命[之]树") == "生命树"
end

@testset "parse_newick" begin
	@test TOL.parse_newick("(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);") ==
		['(', "A", 0.1, "B", 0.2, '(', "C", 0.3, "D", 0.4, ')', 0.5, ')']
	@test TOL.parse_newick("(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;") ==
		['(', "A", 0.1, "B", 0.2, 
			  '(', "C", 0.3, "D", 0.4, ')', ("E",), 0.5, ')', ("F",)]
end

@testset "fromnewick" begin
	global tree = fromnewick("(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;")
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

@testset "gettips" begin
	tree = fromnewick("(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;")
	@test TOL.gettips(tree) == [2, 3, 5, 6]
end
