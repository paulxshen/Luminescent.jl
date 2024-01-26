using AbbreviatedStackTraces, ReverseStackTraces, Pkg
# ENV["JULIA_PKG_PRECOMPILE_AUTO"]
# ENV["JULIA_PKG_PRECOMPILE_AUTO"] = 0
Pkg.UPDATED_REGISTRY_THIS_SESSION[] = true
# using 

Base.round(x::AbstractFloat) = Base.round(Int, x)

T = Float32

Base.getindex(s::Symbol, i) = Symbol(String(s)[i])


f32(x::Real) = Float32(x)
f32(x::AbstractArray) = f32.(x)
° = π / 180