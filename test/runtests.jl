# test/runtests.jl

using TreeOfLife
const TOL = TreeOfLife
using Test

import Aqua

Aqua.test_all(TreeOfLife)

include("types.jl")
include("newick.jl")
include("tree.jl")
