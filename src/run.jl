include("main.jl")
gfrun(lastrun())
# # # using Pkg
# # pkg"add Porcupine,Jello,ArrayPadding;up"

# using Pkg
# pkg"dev C:\Users\pxshe\OneDrive\Desktop\Porcupine.jl;dev C:\Users\pxshe\OneDrive\Desktop\ArrayPadding.jl; dev C:\Users\pxshe\OneDrive\Desktop\Jello.jl;up"
# using FFTW
# a = ones(100, 100)
# @time diff(a, dims=1)
# @time fft(a)
# 1