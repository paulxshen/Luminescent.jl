using Documenter
include("../src/FDTDEngine.jl")
using .FDTDEngine

makedocs(
    sitename="FDTDEngine.jl",
    format=Documenter.HTML(),
    # modules=[FDTDEngine],
    pages=[
        "index.md",
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
run(`cmd /C MOVE /Y makedocs/build/\* docs/`)
run(`cmd /C MOVE /Y examples/\*/\*.mp4 docs/assets/`)