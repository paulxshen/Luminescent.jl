ENV["JULIA_PKG_PRECOMPILE_AUTO"] = 0
include("main.jl")
Base.convert(::Type{Float64}, x::ComplexF64) = real(x)

# picrun(lastrun())#; eta=0.01, iters=20)
# gfrun(lastrun(); gpu_backend="CUDA")
# picrun(lastrun(name="tiny"))
# picrun(lastrun(name="straight"); N=3)
# picrun(lastrun(name="bend"); N=2)
picrun(lastrun(name="bend"); N=3)
# picrun(lastrun(name="1x2_splitter"))
# # # using Pkg
# # pkg"add Porcupine,Jello,ArrayPadding;up"
