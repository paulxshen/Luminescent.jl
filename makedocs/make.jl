using Documenter
include("../src/fdtd_prerelease.jl")
using .fdtd_prerelease

makedocs(
    sitename="fdtd_prerelease",
    format=Documenter.HTML(),
    # modules=[fdtd_prerelease],
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
run(`mv makedocs/build docs`)