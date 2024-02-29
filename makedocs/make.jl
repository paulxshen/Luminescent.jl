using Documenter
include("../src/FDTDEngine.jl")
using .FDTDEngine

makedocs(
    sitename="FDTDEngine.jl",
    format=Documenter.HTML(),
    # modules=[FDTDEngine],
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