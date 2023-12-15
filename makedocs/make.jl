using Documenter
using _fdtd

makedocs(
    sitename = "_fdtd",
    format = Documenter.HTML(),
    modules = [_fdtd]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
