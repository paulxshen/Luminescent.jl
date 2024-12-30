# using Luminescent
include("main.jl")
ENV["JULIA_SSL_CA_ROOTS_PATH"] = ""
ENV["JULIA_PKG_PRECOMPILE_AUTO"] = 0

# using CUDA
# picrun(joinpath("runs", "straight");)# gpuarray=cu)
# picrun(joinpath("runs", "bend_R5"))
# picrun(joinpath("runs", "mode_converter"))
picrun(joinpath("runs", "demux"))
# picrun(joinpath("runs", "splitter"))

