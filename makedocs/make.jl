using Documenter
include("../src/Luminescent.jl")
using .Luminescent

makedocs(
    sitename="Luminescent.jl",
    format=Documenter.HTML(),
    # modules=[Luminescent],
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