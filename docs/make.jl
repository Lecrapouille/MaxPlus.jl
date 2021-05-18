using Documenter, MaxPlus

makedocs(
    modules = [MaxPlus],
    format = Documenter.HTML(),
    pages = Any[
 #       "Home" => "../README.md",
        "Introduction" => "functions.md",
        "Bibliography" => "bibliography.md",
    ],
    sitename = "MaxPlus.jl",
    repo="https://github.com/Lecrapouille/MaxPlus.jl/blob/{commit}{path}#L{line}",
    authors="Quentin Quadrat [AT gmail. COM]"
)

deploydocs(
    repo   = "github.com/Lecrapouille/MaxPlus.jl.git",
)

