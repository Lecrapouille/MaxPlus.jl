using Documenter, SparseArrays, MaxPlus

push!(LOAD_PATH, "../src/")

# We are in the folder MaxPlus.jl/docs/src. Copy files such as MaxPlus.jl/tutorial/README.md
# in this folder else Documenter.jl does not find them.
cp(normpath(@__FILE__, "../../tutorial/README.md"), normpath(@__FILE__, "../src/tutorial.md"); force=true)

# Replace local markdown links to Documenter local links
infile = normpath(@__FILE__, "../../README.md")
outfile = normpath(@__FILE__, "../src/index.md")
out = open(outfile, "w+")
for line in readlines(infile)
    line = replace(line, "](tutorial)" => "](tutorial.md)")
    line = replace(line, "docs/src/" => "")
    write(out, line * "\n")
end
close(out)

# Call Documenter.jl
makedocs(
    modules = [MaxPlus],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        edit_link = "master"
    ),
    sitename = "MaxPlus.jl",
    repo="https://github.com/Lecrapouille/MaxPlus.jl",
    authors="Quentin Quadrat [AT gmail. COM]",
    pages = Any[
        "Home" => "index.md",
        "API: (max,+) Algebra" => "maxplus.md",
        "API: (min,+) Algebra" => "minplus.md",
        "API: (max,+) Linear system" => "syslin.md",
        "Portage: ScicosLab to Julia" => "portage.md",
        "MaxPlus.jl Tutorials" => "tutorial.md",
        "Non Regression Tests" => "tests.md",
        "Bibliography" => "bibliography.md",
    ],
)

deploydocs(
    repo = "github.com/Lecrapouille/MaxPlus.jl.git",
    branch = "gh-pages",
    devbranch = "master",
    devurl = "master",
)

# Remove copied files
rm(normpath(@__FILE__, "../src/index.md"))
rm(normpath(@__FILE__, "../src/tutorial.md"))
