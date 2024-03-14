° = π / 180
# Base.map(f, c::AbstractSet) = f.(c)
Base.round(x::AbstractFloat) = Base.round(Int, x)
Base.ndims(a) = length(size(a))
# Base.size(a) = size(a.data)
# Base.:÷(x, y) = round(Base.div(x, y))

Base.getindex(s::Symbol, i) = Symbol(String(s)[i])

index(v, dx) = round.(Int, v ./ dx .+ 1)
reindex(i, ratio) = round.(Int, (i .- 1) .* ratio .+ 1)
Base.size(x) = (length(x),)

T = Union{Tuple,AbstractArray,Number}
Base.:-(x::T, y::T) = x .- y
Base.:+(x::T, y::T) = x .+ y
Base.:*(x::Number, y::T) = x .* y
Base.:*(x::T, y::Number) = x .* y
Base.:/(x::Number, y::T) = x ./ y
Base.:/(x::T, y::Number) = x ./ y

# U = Union{AbstractDict}
Base.:+(x::AbstractDict, y) = MyDict([k => (x[k] + y) for (k, y) = zip(_keys(x), (y))])
Base.:-(x::AbstractDict, y) = MyDict([k => (x[k] - y) for (k, y) = zip(_keys(x), (y))])

Base.:+(x::NamedTuple, y) = [(x + y) for (k, x, y) = zip(keys(x), x, values(y))]
Base.:-(x::NamedTuple, y) = [(x - y) for (k, x, y) = zip(keys(x), x, values(y))]
# Base.:+(x::NamedTuple, y) = NamedTuple([k => (x + y) for (k, x, y) = zip(keys(x), x, values(y))])
# Base.:-(x::NamedTuple, y) = NamedTuple([k => (x - y) for (k, x, y) = zip(keys(x), x, values(y))])

d2(x) = round.(x, sigdigits=2)

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

function apply!(p, kw)
    [apply!(p[k], kw[k]) for k = _keys(kw)]
    # [apply(p[k], v) for (k, v) = pairs(kw)]
end
function apply!(p; kw...)
    [apply!(p[k], kw[k]) for k = _keys(kw)]
    # [apply(p[k], v) for (k, v) = pairs(kw)]
end
function apply(p, kw)
    # merge((;), [k => apply(p[k], kw[k]) for k = _keys(kw)])
    MyDict([k => apply(p[k], kw[k]) for k = _keys(kw)])
    # (; [k => apply(p[k], kw[k]) for k = _keys(kw)]...)
    # NamedTuple([k => apply(p[k], kw[k]) for k = _keys(kw)])
    # [apply(p[k], kw[k]) for k = _keys(kw)]
    # [apply(p[k], v) for (k, v) = pairs(kw)]
end
function apply(p; kw...)
    # merge((;), [k => apply(p[k], kw[k]) for k = _keys(kw)])
    MyDict([k => apply(p[k], kw[k]) for k = _keys(kw)])
    # (; [k => apply(p[k], kw[k]) for k = _keys(kw)]...)
    # NamedTuple([k => apply(p[k], kw[k]) for k = _keys(kw)])
    # [apply(p[k], kw[k]) for k = _keys(kw)]
    # [apply(p[k], v) for (k, v) = pairs(kw)]
end
function mark(p, kw)
    MyDict([k => mark(p[k], kw[k]) for k = _keys(kw)])
end
function mark(p; kw...)
    # [k => mark(p[k], kw[k]) for k = _keys(kw)]
    # [mark(p[k], kw[k]) for k = _keys(kw)]
    MyDict([k => mark(p[k], kw[k]) for k = _keys(kw)])
end
function unmark(kw)
    MyDict([k => Array(kw[k]) for k = _keys(kw)])
end
function unmark(; kw...)
    MyDict([k => Array(kw[k]) for k = _keys(kw)])
    # [Array(kw[k]) for k = _keys(kw)]
    # [k => Array(kw[k]) for k = _keys(kw)]
end
function mark(v::AbstractVector, a)
    l = sum(v) do p
        p.l
    end
    r = sum(v) do p
        p.r
    end
    PaddedArray(a, l, r)
end

function MyDict(v)
    r = OrderedDict{Symbol,Any}()
    # r = Dict{Symbol,Any}()
    for (k, v) = v
        r[k] = v
    end
    r
end
function _keys(x::OrderedDict)
    k = 0
    ignore_derivatives() do
        k = keys(x) |> collect
    end
    k
end
_keys(x) = keys(x)
footer = "Created by Paul Shen pxshen@alumni.stanford.edu"