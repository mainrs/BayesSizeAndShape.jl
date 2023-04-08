# make.jl
#using Pkg
#Pkg.activate("/Users/gianlucamastrantonio/Dropbox (Politecnico di Torino Staff)/lavori/gitrepo/BayesSizeAndShape/todelete/rat")
using Documenter, BayesSizeAndShape

makedocs(
    modules = [BayesSizeAndShape],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true"),
    authors = "gianluca.mastrantonio",
    sitename = "BayesSizeAndShape.jl",
    pages = Any["index.md"]
    # strict = true,
    # clean = true,
    # checkdocs = :exports,
)

deploydocs(
    repo = "github.com/GianlucaMastrantonio/BayesSizeAndShape.jl.git",
    push_preview = true
)
