include("main.jl")
config_logger(printToStdOut=true, maxLogs=3)
Base.Integer(x::TrackedFloat32) = Int(x)
TrackedFloats.tf_untrack_complex(x) = x
TrackedFloats.tf_track_complex(x) = x
Base.zero(x::Type{Any}) = 0.0f0
Base.sincos(x::TrackedFloat32) = sincos(Float32(x))
Base.:/(x::Float32, y::Float32) = iszero(y) ? error("nan") : Float32(Float64(x) / Float64(y))
Base.:^(x::Float32, y::Float32) = iszero(y) || x == 0 ? error("nan") : Float32(Float64(x)^Float64(y))
gfrun(lastrun())
# # # using Pkg
# # pkg"add Porcupine,Jello,ArrayPadding;up"

# using Pkg
# pkg"dev C:\Users\pxshe\OneDrive\Desktop\Porcupine.jl;dev C:\Users\pxshe\OneDrive\Desktop\ArrayPadding.jl; dev C:\Users\pxshe\OneDrive\Desktop\Jello.jl;up"
7768519