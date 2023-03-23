# docs/make.jl

using TreeOfLife

using Documenter

makedocs(
	sitename = "TreeOfLife.jl",
	pages = [
		"TreeOfLife.jl" => "index.md",
		],
	modules = [TreeOfLife]
)

deploydocs(
    repo = "github.com/Mikumikunisiteageru/TreeOfLife.jl.git",
	versions = ["stable" => "v^", "v#.#.#"]
)
