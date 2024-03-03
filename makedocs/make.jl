using Documenter
include("../src/Luminesce.jl")
using .Luminesce

makedocs(
    sitename="Luminesce.jl",
    format=Documenter.HTML(),
    # modules=[Luminesce],
    pages=[
        "index.md",
        "guide.md",
        "tutorials.md",
        "people.md",
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#