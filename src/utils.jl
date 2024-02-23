° = π / 180

Base.round(x::AbstractFloat) = Base.round(Int, x)
Base.ndims(b::Zygote.Buffer) = Base.ndims(b.data)

Base.getindex(s::Symbol, i) = Symbol(String(s)[i])

index(v, dx) = round.(Int, v ./ dx .+ 1)
reindex(i, ratio) = round.(Int, (i .- 1) .* ratio .+ 1)
Base.size(x::Union{NamedTuple,Tuple}) = (length(x),)

T = Union{Tuple,AbstractArray,Number}
Base.:-(x::T, y::T) = x .- y
Base.:+(x::T, y::T) = x .+ y


# __precompile__(false)
(m::Number)(a...) = m
gaussian(x; μ=0, σ=1) = exp(-((x - μ) / σ)^2)

function place(a, b, o; lazy=false)
    a + pad(b, 0, Tuple(o) .- 1, size(a) .- size(b) .- Tuple(o) .+ 1)
    # place!(buf, b, o)
end
# function place!(a::AbstractArray, b, o)
#     buf = bufferfrom(a)
#     buf = place!(buf, b, o)
#     copy(buf)
# end
function place!(a, b, o)
    a[[i:j for (i, j) = zip(o, o .+ size(b) .- 1)]...] = b
    a
end

function place!(a, b; o)
    # @show size(a), size(b), o
    place!(a, b, o .- floor.((size(b) .- 1) .÷ 2))
end