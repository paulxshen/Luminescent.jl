ENV["JULIA_PKG_PRECOMPILE_AUTO"] = 0
include("main.jl")
Base.convert(::Type{Float64}, x::ComplexF64) = real(x)

# picrun(lastrun(); eta=0.4)
# gfrun(lastrun(); gpu_backend="CUDA")
# picrun(lastrun("tiny"))
# picrun(lastrun("straight"); N=3)
# picrun(lastrun("simtest2"; wd=joinpath("build", "precompile_execution")))
# picrun(lastrun("bend"); N=2)
picrun(lastrun("bend"); N=3)
# picrun(lastrun("1x2_splitter"))
# # # using Pkg
# # pkg"add Porcupine,Jello,ArrayPadding;up"
0
# heatmap(abs.(sum(res.source_instances[1].sigmodes[1][2].Ez, dims=1)[1, :, :]))