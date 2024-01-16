using AbbreviatedStackTraces, ReverseStackTraces, Pkg
ENV["JULIA_PKG_PRECOMPILE_AUTO"]
ENV["JULIA_PKG_PRECOMPILE_AUTO"] = 0
Pkg.UPDATED_REGISTRY_THIS_SESSION[] = true
# using 

Base.round(x::AbstractFloat) = Base.round(Int, x)

T = Float32

# Base.eachslice(a) = eachslice(a, dims=ndims(a))

Base.getindex(s::Symbol, i) = Symbol(String(s)[i])

Base.:+(x::AbstractArray, y::Number) = x .+ y
Base.:+(x::Number, y::AbstractArray) = x .+ y
Base.:-(x::AbstractArray, y::Number) = x .- y
Base.:-(x::Number, y::AbstractArray) = x .- y