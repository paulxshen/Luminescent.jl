ENV["JULIA_PKG_PRECOMPILE_AUTO"] = 0
include("main.jl")
# using GLMakie

# picrun(lastrun())#; eta=0.01, iters=20)
# gfrun(lastrun(); gpu_backend="CUDA")
picrun(lastrun(name="1x2_splitter"))
# picrun(lastrun(name="tiny"))
# picrun(lastrun(name="straight");)
# # # using Pkg
# # pkg"add Porcupine,Jello,ArrayPadding;up"

# using Pkg
# pkg"dev C:\Users\pxshe\OneDrive\Desktop\beans\Porcupine.jl;dev C:\Users\pxshe\OneDrive\Desktop\beans\ArrayPadding.jl; dev C:\Users\pxshe\OneDrive\Desktop\beans\Jello.jl;dev C:\Users\pxshe\OneDrive\Desktop\beans\VectorModesolver.jl;up"
# pkg"add https://github.com/hammy4815/VectorModesolver.jl"
# pkg"dev ~/Porcupine.jl;dev ~/ArrayPadding.jl; dev ~/Jello.jl;up"
# 7768519
# using CairoMakie
# heatmap(abs.(mode_solutions[1].calibrated_modes[1].Ex))
# # volume(abs.(sols[1].u.Ey))
# using GLMakie
# volume(abs.(run_probs[1].source_instances[1]._g.Jy))
# heatmap(sum(abs.(run_probs[1].source_instances[1]._g.Jy), dims=1)[1, :, :])

