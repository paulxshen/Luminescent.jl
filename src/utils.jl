function place_mask(a, o, mask, c=1)
    o = round.(o)
    _mask = (mask / maximum(abs, mask)) .> 0.5
    _ksam = .!(_mask)
    a = a .* pad(_ksam, true, o - 1, size(a) .- size(_ksam) - o + 1)
    place(a, o, _mask .* mask,)
end
d2(x) = round.(x, sigdigits=2)

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

# Flux.gpu(d::T) where {T<:Dictlike} = dict(T, [k => Flux.gpu(d[k]) for k = keys(d)])
# Flux.cpu(d::T) where {T<:Dictlike} = dict(T, [k => Flux.cpu(d[k]) for k = keys(d)])
Flux.gpu(d::Dictlike) = apply(Flux.gpu, d)
Flux.cpu(d::Dictlike) = apply(Flux.cpu, d)
Flux.gpu(v::AbstractVector{T}) where {T<:Dictlike} = gpu.(v)
Flux.cpu(v::AbstractVector{T}) where {T<:Dictlike} = cpu.(v)
# Base.getproperty(d::AbstractDict, k::Symbol) = hasfield(d, k) ? getfield(d, k) : (haskey(d, k) ? d[k] : d[string(k)])
Base.getproperty(d::AbstractDict, k::Symbol) = hasproperty(d, k) ? getfield(d, k) : d[k]

footer = "Suppress this message by verbose=false\n2024 (c) Paul Shen at Luminescent AI\n<pxshen@alumni.stanford.edu>"