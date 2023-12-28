using AbbreviatedStackTraces
ENV["JULIA_PKG_PRECOMPILE_AUTO"] = 0
# using ReverseStackTraces

Base.round(x::AbstractFloat) = Base.round(Int, x)

# T = Float32
# Base.:+(x::Float64, y::Float32) = Float32(x) + y
# Base.:+(x::Float32, y::Float64) = x + Float32(y)
# Base.:-(x::Float64, y::Float32) = Float32(x) - y
# Base.:-(x::Float32, y::Float64) = x - Float32(y)
# Base.:*(x::Float64, y::Float32) = Float32(x) * y
# Base.:*(x::Float32, y::Float64) = x * Float32(y)
# Base.:/(x::Float64, y::Float32) = Float32(x) / y
# Base.:/(x::Float32, y::Float64) = x / Float32(y)

Base.eachslice(a) = eachslice(a, dims=ndims(a))

Base.getindex(s::Symbol, i) = Symbol(String(s)[i])