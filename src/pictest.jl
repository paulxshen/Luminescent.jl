# using Luminescent
include("main.jl")
ENV["JULIA_SSL_CA_ROOTS_PATH"] = ""
ENV["JULIA_PKG_PRECOMPILE_AUTO"] = 0
using CUDA
Random.seed!(1234)

# picrun(joinpath("runs", "straight");)# array=cu)
# picrun(joinpath("runs", "bend_R5"), array=cu)
# picrun(joinpath("runs", "mode_converter"))

# picrun(joinpath("runs", "splitter"); array=cu)
# picrun(joinpath("runs", "splitter"))

picrun(joinpath("build", "precompile_execution", "tiny"))
picrun(joinpath("build", "precompile_execution", "tinycu"))
# picrun(joinpath("runs", "tiny3"))
# picrun(joinpath("runs", "back"))# array=cu)
# models[1]()

# picrun(joinpath("runs", "demux"))
# picrun(joinpath("runs", "straight"))