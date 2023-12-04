using AbbreviatedStackTraces
# using ReverseStackTraces

Base.round(x::AbstractFloat) = Base.round(Int, x)

T = Float32
Base.:+(x::Float64, y::Float32) = Float32(x) + y
Base.:+(x::Float32, y::Float64) = x + Float32(y)
Base.:-(x::Float64, y::Float32) = Float32(x) - y
Base.:-(x::Float32, y::Float64) = x - Float32(y)
Base.:*(x::Float64, y::Float32) = Float32(x) * y
Base.:*(x::Float32, y::Float64) = x * Float32(y)
Base.:/(x::Float64, y::Float32) = Float32(x) / y
Base.:/(x::Float32, y::Float64) = x / Float32(y)

Base.getindex(s::Symbol, i) = Symbol(String(s)[i])