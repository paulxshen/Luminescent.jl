# using Luminescent
include("main.jl")
ENV["JULIA_SSL_CA_ROOTS_PATH"] = ""
ENV["JULIA_PKG_PRECOMPILE_AUTO"] = 0
using CUDA


# picrun(joinpath("runs", "straight");)# array=cu)
# picrun(joinpath("runs", "bend_R5"), array=cu)
# picrun(joinpath("runs", "mode_converter"))
# picrun(joinpath("runs", "demux"))

# picrun(joinpath("runs", "splitter"); array=cu)
# picrun(joinpath("runs", "splitter"))

# picrun(joinpath("runs", "back");)# array=cu)
picrun(joinpath("runs", "tiny"); array=cu)
# models[1]()
