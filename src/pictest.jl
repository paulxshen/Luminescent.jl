ENV["JULIA_PKG_PRECOMPILE_AUTO"] = 0
include("main.jl")

picrun(lastrun())#; gpu_backend="CUDA")

# picrun(lastrun("bend_R5"))
# picrun(lastrun("bend_R5_multi"))
# picrun(lastrun("straight"); N=2)
# picrun(lastrun("bend_R2__1.1"); N=2)
# picrun(lastrun("bend_R2_multi"); N=2)
# picrun(lastrun("tiny"))

# picrun(lastrun("mode_converter");)# eta=0.4)
# picrun(lastrun("simtest2"; wd=joinpath("build", "precompile_execution")))
# picrun(lastrun("bend"); N=2)
# picrun(lastrun("splitter"))
# # # using Pkg
# # pkg"add Porcupine,Jello,ArrayPadding;up"
# heatmap(abs.(sum(res.source_instances[1].sigmodes[1][2].Ez, dims=1)[1, :, :]))