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

function apply(d1, d2; kw...)
    dict([k => apply(d1[k], d2[k]; kw...) for k = keys(d2)])
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


# Flux.gpu(d::T) where {T<:Dictlike} = dict(T, [k => Flux.gpu(d[k]) for k = keys(d)])
# Flux.cpu(d::T) where {T<:Dictlike} = dict(T, [k => Flux.cpu(d[k]) for k = keys(d)])

# Base.Float16(x) = f16(x)
# Base.Float32(x) = f32(x)


footer = "Suppress this message by verbose=false\n2024 (c) Paul Shen at Luminescent AI\n<pxshen@alumni.stanford.edu>"