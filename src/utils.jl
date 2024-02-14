° = π / 180

Base.round(x::AbstractFloat) = Base.round(Int, x)
Base.ndims(b::Zygote.Buffer) = Base.ndims(b.data)

Base.getindex(s::Symbol, i) = Symbol(String(s)[i])

index(v, dx) = round.(Int, v ./ dx .+ 1)
reindex(i, ratio) = round.(Int, (i .- 1) .* ratio .+ 1)
Base.size(x::NamedTuple) = (length(x),)

T = Union{Tuple,AbstractArray,Number}
Base.:-(x::T, y::T) = x .- y
Base.:+(x::T, y::T) = x .+ y

Base.view(b::Buffer, i...) = b[i...]

f32(x::Real) = Float32(x)
f32(x::AbstractArray) = f32.(x)

# __precompile__(false)
(m::Number)(a...) = m
gaussian(x; μ=0, σ=1) = exp(-((x - μ) / σ)^2)

function place(a, b, start; lazy=false)
    a + pad(b, 0, Tuple(start) .- 1, size(a) .- size(b) .- Tuple(start) .+ 1)
    # place!(buf, b, start)
end
function place!(a::AbstractArray, b, start)
    buf = bufferfrom(a)
    place!(buf, b, start)
    copy(buf)
end
function place!(a, b, start)
    a[[i:j for (i, j) = zip(start, start .+ size(b) .- 1)]...] = b
end

function place!(a, b; o)
    # @show size(a), size(b), o
    place!(a, b, o .- floor.((size(b) .- 1) .÷ 2))
end