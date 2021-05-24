using Documenter, SparseArrays, MaxPlus

# Hack Documenter.jl seems not finding files outside docs/ folder: so copy them
cp(normpath(@__FILE__, "../../README.md"), normpath(@__FILE__, "../src/home.md"); force=true)
cp(normpath(@__FILE__, "../../tutorial/README.md"), normpath(@__FILE__, "../src/tutorials.md"); force=true)

makedocs(
    modules = [MaxPlus],
    format = Documenter.HTML(),
    sitename = "MaxPlus.jl",
#   repo="https://github.com/Lecrapouille/MaxPlus.jl/blob/{commit}{path}#L{line}",
    repo="https://github.com/Lecrapouille/MaxPlus.jl",
    authors="Quentin Quadrat [AT gmail. COM]",
    pages = Any[
        "Home" => "home.md",
        "API" => "API.md",
        "Tutorials" => "tutorials.md",
        "ScicosLab to Julia" => "functions.md",
        "Bibliography" => "bibliography.md",
    ],
)

deploydocs(
    repo = "github.com/Lecrapouille/MaxPlus.jl.git",
    branch = "gh-pages",
    devbranch = "dev",
    devurl = "dev",
)

# Remove copied files
rm(normpath(@__FILE__, "../src/home.md"))
rm(normpath(@__FILE__, "../src/tutorials.md"))
