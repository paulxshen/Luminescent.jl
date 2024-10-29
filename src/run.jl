include("main.jl")
# gfrun(lastrun(); gpu_backend="CUDA")
picrun(lastrun(name="bend"))
# picrun(lastrun(name="straight");)
# # # using Pkg
# # pkg"add Porcupine,Jello,ArrayPadding;up"

# using Pkg
# pkg"dev C:\Users\pxshe\OneDrive\Desktop\Porcupine.jl;dev C:\Users\pxshe\OneDrive\Desktop\ArrayPadding.jl; dev C:\Users\pxshe\OneDrive\Desktop\Jello.jl;up"
# 7768519
# using CairoMakie
# heatmap(abs.(mode_solutions[1].calibrated_modes[1].Ex))
# # volume(abs.(sols[1].u.Ey))
# using GLMakie
# volume(abs.(run_probs[1].source_instances[1]._g.Jy))
# heatmap(sum(abs.(run_probs[1].source_instances[1]._g.Jy), dims=1)[1, :, :])
# heatmap(abs.(_sources[1].mode.