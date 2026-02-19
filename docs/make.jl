using Documenter
using NURBS

makedocs(
    sitename = "NURBS.jl",
    modules = [NURBS],
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
