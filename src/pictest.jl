# using Luminescent
include("main.jl")
ENV["JULIA_SSL_CA_ROOTS_PATH"] = ""
ENV["JULIA_PKG_PRECOMPILE_AUTO"] = 0
Random.seed!(1234)

# picrun(joinpath("runs", "straight");)# array=cu)
# picrun(joinpath("runs", "bend_R5"), array=cu)
# picrun(joinpath("runs", "mode_converter"))

# picrun(joinpath("runs", "splitter"); array=cu)
picrun(joinpath("runs", "splitter"))
# using GLMakie: volume
# picrun(joinpath("build", "precompile_execution", "tiny_2_float32_CUDA"))#; framerate=10)
# picrun(joinpath("build", "precompile_execution", "tiny_3_float32_CUDA"))
# picrun(joinpath("build", "precompile_execution", "tiny_3_float32_None"))
# picrun(joinpath("build", "precompile_execution", "back_float32"))
# picrun(joinpath("runs", "tiny3"))
# picrun(joinpath("runs", "back"))# array=cu)
# models[1]()

# picrun(joinpath("runs", "demux"))
# picrun(joinpath("runs", "straight"))

# for p = readdir("build/precompile_execution", join=true)
#     # if !contains(string(p), "16") && !contains(string(p), "back")
#     if !contains(string(p), "16") #&& contains(string(p), "back")
#         # if !contains(string(p), "16") && contains(string(p), "back")
#         picrun(p)
#     end
# end
# # 
