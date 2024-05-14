° = π / 180
gaussian(x; μ=0, σ=1) = exp(-((x - μ) / σ)^2)

d2(x) = round.(x, sigdigits=2)

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

function place!(a, b; center)
    # @show size(a), size(b), o
    place!(a, b, center .- floor.((size(b) .- 1) .÷ 2))
end

function apply!(p, kw)
    dict([k => apply!(p[k], kw[k]) for k = keys(kw)])
    # [apply(p[k], v) for (k, v) = pairs(kw)]
end
function apply!(p; kw...)
    apply!(p, kw)
end
# apply!(p, x::Number) = x
apply(p, x::Number) = apply(p, x * ones(typeof(x), p.default_size))

function apply(p, kw)
    dict([k => apply(p[k], kw[k]) for k = keys(kw)])
    # [apply(p[k], v) for (k, v) = pairs(kw)]
end
function apply(p; kw...)
    apply(p, kw)
end

function apply(d::Dictlike, a::AbstractArray)
    dict([k => apply(d[k], a) for k = keys(d)])
end

function apply(i::AbstractVector{<:AbstractRange}, a::AbstractArray)
    getindex(a, i...)
end
function mark(p, kw)
    dict([k => mark(p[k], kw[k]) for k = keys(kw)])
end
function mark(p; kw...)
    # [k => mark(p[k], kw[k]) for k = keys(kw)]
    # [mark(p[k], kw[k]) for k = keys(kw)]
    dict([k => mark(p[k], kw[k]) for k = keys(kw)])
end
function unmark(kw)
    dict([k => Array(kw[k]) for k = keys(kw)])
end
function unmark(; kw...)
    dict([k => Array(kw[k]) for k = keys(kw)])
    # [Array(kw[k]) for k = keys(kw)]
    # [k => Array(kw[k]) for k = keys(kw)]
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

Flux.gpu(d::T) where {T<:Dictlike} = dict(T, [k => Flux.gpu(d[k]) for k = keys(d)])
Flux.cpu(d::T) where {T<:Dictlike} = dict(T, [k => Flux.cpu(d[k]) for k = keys(d)])
Flux.gpu(v::AbstractVector{T}) where {T<:Dictlike} = gpu.(v)
Flux.cpu(v::AbstractVector{T}) where {T<:Dictlike} = cpu.(v)
# Base.getproperty(d::AbstractDict, k::Symbol) = hasfield(d, k) ? getfield(d, k) : (haskey(d, k) ? d[k] : d[string(k)])
Base.getproperty(d::AbstractDict, k::Symbol) = hasproperty(d, k) ? getfield(d, k) : d[k]

footer = "Suppress this message by verbose=false\n2024 (c) Paul Shen at Luminescent AI\n<pxshen@alumni.stanford.edu>"