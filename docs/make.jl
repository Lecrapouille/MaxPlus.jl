using Pkg

# Get the directory where this script is located
const DOC_DIR = dirname(@__FILE__)

# Activate and instantiate docs environment
Pkg.activate(DOC_DIR)
Pkg.instantiate()

# Set GR to no-display mode for CI/build documentation
ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")

# Import necessary packages
using Documenter
using SparseArrays
using MaxPlus

# Make sure local src directory is in the load path
push!(LOAD_PATH, joinpath(DOC_DIR, "..", "src"))

# Copy the tutorial's README.md for Documenter's use in the docs
cp(joinpath(DOC_DIR, "..", "tutorial", "README.md"),
   joinpath(DOC_DIR, "src", "tutorial.md"); force=true)

# Generate index.md from top-level README.md, adjusting links for Documenter
infile = joinpath(DOC_DIR, "..", "README.md")
outfile = joinpath(DOC_DIR, "src", "index.md")
open(outfile, "w") do out
    for line in readlines(infile)
        # Fix tutorial links and strip "docs/src/" prefix from local references
        line = replace(line, "](tutorial)" => "](tutorial.md)")
        line = replace(line, "docs/src/" => "")
        write(out, line * "\n")
    end
end

makedocs(
    modules = [MaxPlus],
    # Only warn on @docs/reference errors (API pages may be out of sync)
    warnonly = true,
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        edit_link = "master",
    ),
    sitename = "MaxPlus.jl",
    repo = "https://github.com/Lecrapouille/MaxPlus.jl",
    authors = "Quentin Quadrat [AT gmail. COM]",
    pages = [
        "Home" => "index.md",
        "API: (max,+) Algebra" => "maxplus.md",
        "API: (min,+) Algebra" => "minplus.md",
        "API: (max,+) Linear system" => "syslin.md",
        "Portage: ScicosLab to Julia" => "portage.md",
        "Flowshop Example" => "flowshop.md",
        "MaxPlus.jl Tutorials" => "tutorial.md",
        "Non Regression Tests" => "tests.md",
        "Bibliography" => "bibliography.md",
    ]
)

deploydocs(
    repo = "github.com/Lecrapouille/MaxPlus.jl.git",
    branch = "gh-pages",
    devbranch = "master",
    devurl = "master",
)

# Clean up generated files
rm(joinpath(DOC_DIR, "src", "index.md"); force=true)
rm(joinpath(DOC_DIR, "src", "tutorial.md"); force=true)
