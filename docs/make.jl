using Documenter, SparseArrays, MaxPlus

push!(LOAD_PATH, "../src/")

# We are in the folder MaxPlus.jl/docs/src. Copy files such as MaxPlus.jl/tutorial/README.md
# in this folder else Documenter.jl does not find them.
cp(normpath(@__FILE__, "../../README.md"), normpath(@__FILE__, "../src/index.md"); force=true)
cp(normpath(@__FILE__, "../../tutorial/README.md"), normpath(@__FILE__, "../src/tutorial.md"); force=true)

makedocs(
    modules = [MaxPlus],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        edit_link = "dev"
    ),
    sitename = "MaxPlus.jl",
    repo="https://github.com/Lecrapouille/MaxPlus.jl",
    authors="Quentin Quadrat [AT gmail. COM]",
    pages = Any[
        "Home" => "index.md",
        "API" => "API.md",
        "Tutorials" => "tutorial.md",
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
rm(normpath(@__FILE__, "../src/index.md"))
rm(normpath(@__FILE__, "../src/tutorial.md"))
