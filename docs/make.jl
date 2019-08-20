if VERSION < v"1.0.1"
    # Workaround for JuliaLang/julia/pull/28625
    if Base.HOME_PROJECT[] !== nothing
        Base.HOME_PROJECT[] = abspath(Base.HOME_PROJECT[])
    end
end

using Documenter, MaxPlus

makedocs(
    modules = [MaxPlus],
    sitename = "MaxPlus",
    format = Documenter.HTML(),
    pages = [
        "Introduction" => "index.md",
    ]
)

deploydocs(
    repo   = "github.com/Lecrapouille/MaxPlus.jl.git",
)

