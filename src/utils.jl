Base.round(x::AbstractFloat) = Base.round(Int, x)
# Base.eachslice(a) = eachslice(a, dims=ndims(a))

Base.getindex(s::Symbol, i) = Symbol(String(s)[i])

index(v, dx) = round.(Int, v ./ dx .+ 1)
reindex(i, ratio) = round.(Int, (i .- 1) .* ratio .+ 1)
Base.size(x::NamedTuple) = (length(x),)

Base.:+(x::AbstractArray, y::Number) = x .+ y
Base.:+(x::Number, y::AbstractArray) = x .+ y
Base.:-(x::AbstractArray, y::Number) = x .- y
Base.:-(x::Number, y::AbstractArray) = x .- y