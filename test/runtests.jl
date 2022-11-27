using TreeOfLife
using Test

@test isa(read_newick("(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);"), Chronogram)
