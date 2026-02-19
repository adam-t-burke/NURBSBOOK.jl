using Documenter
using NURBSBOOK

makedocs(
    sitename = "NURBSBOOK.jl",
    modules = [NURBSBOOK],
    format = Documenter.HTML(prettyurls = false),
    pages = [
        "Home" => "index.md",
        "Types" => "types.md",
        "Basis Functions" => "basis.md",
        "Curves" => "curves.md",
        "Surfaces" => "surfaces.md",
        "Algorithms" => "algorithms.md",
        "API Reference" => "api.md",
    ],
)
