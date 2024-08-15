function resize(a, sz)
    size(a) != Tuple(sz) && @warn "array size $(size(a)) not same as new size $sz. array will be interpolated"
    imresize(a, sz, method=ImageTransformations.Lanczos4OpenCV())

end
"""
"""
struct ModalSource
    f
    center
    lb
    ub
    normal
    tangent
    # zaxis
    # xaxis
    mode
    label
    meta
    function ModalSource(f, mode, center::Base.AbstractVecOrTuple, normal, tangent, lb, ub, ; label::String="", meta=Dict())
        new(f, center, lb, ub, normal, tangent, E2J(mode), label, meta)
    end
    function ModalSource(f, mode, center::Base.AbstractVecOrTuple, normal, tangent, L, ; label::String="", meta=Dict())
        new(f, center, -L / 2, L / 2, normal, tangent, E2J(mode), label, meta)
    end

end
"""
    function PlaneWave(f, dims; fields...)

Constructs plane wave source

Args
- f: time function
- dims: eg -1 for wave coming from -x face
- fields: which fields to excite & their scaling constants (typically a current source, eg Jz=1)
"""
struct PlaneWave
    f
    fields
    dims
    label
    function PlaneWave(f, dims, label=""; fields...)
        new(f, fields, dims, label)
    end
end
@functor PlaneWave
"""
    function GaussianBeam(f, σ, c, dims; fields...)

Constructs gaussian beam source

Args
- f: time function
- σ: std dev length
- dims: eg 1 for x direction
- fields: which fields to excite & their scaling constants (typically a current source, eg Jz=1)
"""
struct GaussianBeam
    f
    σ
    fields
    c
    dims
    function GaussianBeam(f, σ, c, dims; fields...)
        new(f, σ, fields, c, dims)
    end
end
@functor GaussianBeam

"""
    function Source(f, c, lb, ub, label=""; fields...)
    function Source(f, c, L, label=""; fields...)
        
Constructs custom  source. Can be used to specify uniform or modal sources

Args
- f: time function
- c: origin or center of source
- lb: lower bounds wrt to c
- ub: upper bounds wrt to c
- L: source dimensions in [wavelengths]
- fields: which fields to excite & their scaling constants (typically a current source, eg Jz=1)
"""
struct Source
    f
    fields
    center
    lb
    ub
    label
    meta
    function Source(f, fields, c, lb, ub, label::String=""; kw...)
        new(f, fields, c, lb, ub, label, kw)
    end
    function Source(f, fields, c, L, label::String=""; kw...)
        new(f, fields, c, -L / 2, L / 2, label, kw)
    end
end
# fields(m::Source) = m.fields
Base.string(m::Union{Source,ModalSource}) =
    """
    $(m.label): $(count((m.ub.-m.lb).!=0))-dimensional source, centered at $(m.center|>d2), spanning from $(m.lb|>d2) to $(m.ub|>d2) relative to center,"""
#  exciting $(join(keys(m.fields),", "))"""

function Source(f, c, L, label::AbstractString=""; fields...)
    # Source
    # function Source(f, c, L::Union{AbstractVector{<:Real},Tuple{<:Real}}; fields...)
    Source(f, c, -L / 2, L / 2, label; fields...)
end
Source = Source

struct SourceInstance
    f
    k
    g
    _g
    o
    center
    label
    meta
end
@functor SourceInstance (g, _g,)

function SourceInstance(s::ModalSource, dx, sizes, common_left_pad_amount, fl, sz0; F=Float32)
    @unpack f, center, lb, ub, normal, tangent, meta = s
    C = complex(F)
    J = dict([:Jx => [C(0)], :Jy => [C(0)], :Jz => [C(0)]])
    for k = keys(s.mode)
        J[k] = s.mode[k]
    end
    J = values(J)

    L = ub .- lb
    sz = max.(1, round.(Int, L ./ dx)) |> Tuple
    d = length(lb)
    D = length(center) # 2D or 3D
    if D == 2
        zaxis = [normal..., 0]
        yaxis = [0, 0, 1]
        # xaxis = cross(yaxis, zaxis)
        xaxis = cross(yaxis, zaxis)
        # @show J
        # J = reshape.([Jx, Jy, Jz], (sz,))
    else
        zaxis = normal |> F
        xaxis = tangent |> F
        yaxis = cross(zaxis, xaxis)
    end

    frame = [xaxis, yaxis, zaxis]
    # J = resize.(J, (sz,))
    J = reframe(frame, J)
    if D == 2
        J = [J[:, :, 1] for J in J]
    end
    lb = sum(lb .* frame[1:d])[1:D]
    ub = sum(ub .* frame[1:d])[1:D]
    L = ub - lb

    Jx, Jy, Jz = J
    if D == 2
        fields = (; Jx, Jy)
    else
        fields = (; Jx, Jy, Jz)
    end
    v = zip(lb, ub)
    lb = minimum.(v)
    ub = maximum.(v)
    SourceInstance(Source(f, fields / dx^(D - d), center, lb, ub; meta...), dx, sizes, common_left_pad_amount, fl, sz0; F)
end
function SourceInstance(s::PlaneWave, dx, sizes, common_left_pad_amount, fl, sz0; F=Float32)
    @unpack f, fields, dims, label, meta = s
    _F(x::Real) = F(x)
    _F(x::Complex) = complex(F)(x)
    f = _F ∘ f
    d = length(common_left_pad_amount)
    g = Dict([k => _F(fields[k]) * ones([i == abs(dims) ? 1 : sz0[i] for i = 1:d]...) / dx for k = keys(fields)])
    o = NamedTuple([k =>
        1 .+ fl[k] .+ (dims < 0 ? 0 : [i == abs(dims) ? sizes[k][i] - 1 : 0 for i = 1:d])
                    for k = keys(fl)])
    _g = Dict([k => place(zeros(F, sizes[k]), o[k], g[k],) for k = keys(fields)])
    c = first(values(sizes)) .÷ 2
    SourceInstance(f, keys(fields), g, _g, o, c, label, meta)
end

function SourceInstance(s::GaussianBeam, dx, sizes, fl, stop; F=Float32)
    _F(x::Real) = F(x)
    _F(x::Complex) = complex(F)(x)
    f = _F ∘ f
    @unpack f, σ, fields, c, dims = s
    n = round(Int, 2σ / dx)
    r = n * dx
    r = [i == abs(dims) ? (0:0) : range(-r, r, length=(2n + 1)) for i = 1:length(c)]
    g = [gaussian(norm(F.(collect(v)))) for v = Iterators.product(r...)] / dx
    fl = fl .- 1 .+ index(c, dx) .- round.(Int, (size(g) .- 1) ./ 2)
    _g = place(zeros(F, sz), g, fl)
    SourceInstance(f, keys(fields), g, _g, fl, c, label)
end

function SourceInstance(s::Source, dx, sizes, field_origin, common_left_pad_amount, stop; F=Float32)
    @unpack f, fields, center, lb, ub, label, meta = s
    _F(x::Real) = F(x)
    _F(x::Complex) = complex(F)(x)

    f = _F ∘ f
    g = Dict([k =>
        begin
            if isa(fields[k], AbstractArray)
                # imresize(fields[k], ratio=1)
                sz0 = size(fields[k])
                sz = max.(1, round(abs.(ub - lb) ./ dx)) |> Tuple
                # sz0 != sz && @warn "source array size$sz0 not same  as domain size $sz. source will be interpolated"
                imresize(_F.(fields[k]), sz, method=ImageTransformations.Lanczos4OpenCV())
            else
                r = [a == b ? (a:a) : (a+dx/2*sign(b - a):dx*sign(b - a):b) for (a, b) = zip(lb, ub)]
                [_F.(fields[k](v...)) for v = Iterators.product(r...)]
            end
        end

              for k = keys(fields)])
    o = NamedTuple([k => F((center + lb - o) / dx .+ 1.5) for (k, o) = pairs(field_origin)])
    _g = Dict([k => place(zeros(F, sizes[k]), o[k], g[k],) for k = keys(fields)])
    _center = round(center / dx) + 1 + common_left_pad_amount
    SourceInstance(f, keys(fields), g, _g, o, _center, label, meta)
end

# Complex
apply(v::AbstractVector{<:SourceInstance}, t::Real; kw...) = apply(v, t, kw)
function apply(v::AbstractVector{<:SourceInstance}, t::Real, kw)
    dict([
        # sum([real(s.f(t) .* s._g[k]) for s = v if k in s.k], init=kw[k])
        k => begin

            a = kw[k]
            # gpu = isa(first(v).g, CuArray)
            # if gpu
            #     ignore() do
            #         a = cpu(a)
            #     end
            # end
            for s = v
                if k in s.k
                    mode = s._g[k]
                    # if gpu
                    #     ignore() do
                    #         mode = cpu(mode)
                    #     end
                    # end
                    a += real(s.f(t) .* mode)
                end
            end
            # if gpu
            #     ignore() do
            #         a = gpu(a)
            #     end
            # end
            a
        end
        for k = keys(kw)
        # end for (k, a) = pairs(kw)
    ])
end

# function apply(d,t; kw...)
#     [apply(d[k], kw[k]) for k = keys(kw)]
# end

function E2J(d)
    k0 = filter(k -> string(k)[1] == 'E', keys(d))
    k = [Symbol("J" * string(k)[2:end]) for k in k0]
    v = [d[k] for k in k0]
    dict(Pair.(k, v))
end
function EH2JM(d::T) where {T}
    dict(Pair.(replace(keys(d), :Ex => :Jx, :Ey => :Jy, :Ez => :Jz, :Hx => :Mx, :Hy => :My, :Hz => :Mz), values(d)))
end
# function EH2JM(d::NamedTuple)
#     NamedTuple((d) |> pairs |> OrderedDict |> EH2JM |> pairs)
# end