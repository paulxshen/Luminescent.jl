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
    # zaxis
    # xaxis
    fields
    label
    function ModalSource(f, mode, center::Base.AbstractVecOrTuple, lb, ub, normal, label::String="";)
        new(f, center, lb, ub, normal, E2J(mode), label)
    end
    function ModalSource(f, mode, center::Base.AbstractVecOrTuple, L, normal, label::String="";)
        new(f, center, -L / 2, L / 2, normal, E2J(mode), label)
    end

    function ModalSource(f, center::Base.AbstractVecOrTuple, lb, ub, normal, label::String=""; fields...)
        new(f, center, lb, ub, normal, fields, label)
    end
    function ModalSource(f, center::Base.AbstractVecOrTuple, L, normal, label::String=""; fields...)
        new(f, center, -L / 2, L / 2, normal, fields, label)
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
    function Source(f, c, lb, ub, label::String=""; fields...)
        new(f, fields, c, lb, ub, label)
    end
    function Source(f, c, L, label::String=""; fields...)
        new(f, fields, c, -L / 2, L / 2, label)
    end
end
# fields(m::Source) = m.fields
Base.string(m::Union{Source,ModalSource}) =
    """
    $(m.label): $(count((m.ub.-m.lb).!=0))-dimensional source, centered at $(m.center|>d2), spanning from $(m.lb|>d2) to $(m.ub|>d2) relative to center, exciting $(join(keys(m.fields),", "))"""
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
    c
    label
end
@functor SourceInstance (g, _g,)

function SourceInstance(s::ModalSource, dx, sizes, common_left_pad_amount, fl, sz0; F=Float32)
    @unpack f, center, lb, ub, normal = s
    fields = DefaultDict([0],)
    for k = keys(s.fields)
        fields[k] = s.fields[k]
    end
    @unpack Jx, Jy, Jz = fields
    L = ub .- lb
    sz = max.(1, round.(Int, L ./ dx)) |> Tuple

    d = length(center) # 2D or 3D
    if d == 2
        zaxis = [normal..., 0]
        yaxis = [0, 0, 1]
        # xaxis = cross(yaxis, zaxis)
        xaxis = cross(yaxis, zaxis)
        # @show J
        # J = reshape.([Jx, Jy, Jz], (sz,))
    else
        yaxis = cross(zaxis, xaxis)
    end

    frame = [xaxis, yaxis, zaxis]
    J = [Jx, Jy, Jz]
    # J = resize.(J, (sz,))
    J = reframe(frame, J)
    if d == 2
        J = [J[:, :, 1] for J in J]
    end

    Jx, Jy, Jz = J
    if d == 2
        fields = (; Jx, Jy)
    else
        fields = (; Jx, Jy, Jz)
    end
    SourceInstance(Source(f, center, L; fields...), dx, sizes, common_left_pad_amount, fl, sz0; F)
end
function SourceInstance(s::PlaneWave, dx, sizes, common_left_pad_amount, fl, sz0; F=Float32)
    @unpack f, fields, dims, label = s
    _F(x::Real) = F(x)
    _F(x::Complex) = complex(F)(x)
    f = _F ∘ f
    d = length(common_left_pad_amount)
    g = Dict([k => _F(fields[k]) * ones([i == abs(dims) ? 1 : sz0[i] for i = 1:d]...) / dx for k = keys(fields)])
    o = NamedTuple([k =>
        1 .+ fl[k] .+ (dims < 0 ? 0 : [i == abs(dims) ? sizes[k][i] - 1 : 0 for i = 1:d])
                    for k = keys(fl)])
    _g = Dict([k => place(zeros(F, sizes[k]), g[k], o[k]) for k = keys(fields)])
    c = first(values(sizes)) .÷ 2
    SourceInstance(f, keys(fields), g, _g, o, c, label)
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
    @unpack f, fields, center, lb, ub, label = s
    _F(x::Real) = F(x)
    _F(x::Complex) = complex(F)(x)

    f = _F ∘ f
    r = [a == b ? (a:a) : a+dx/2:dx:b for (a, b) = zip(lb, ub)]
    C = F(1 / dx^count(lb .== ub))
    g = Dict([k =>
        C * begin
            if isa(fields[k], AbstractArray)
                # imresize(fields[k], ratio=1)
                sz0 = size(fields[k])
                sz = Tuple(length.(r))
                sz0 != sz && @warn "source spatial profile array $sz0 not same size as source domain $sz. profile will be interpolated"
                imresize(_F.(fields[k]), sz, method=ImageTransformations.Lanczos4OpenCV())
            else
                [_F.(fields[k](v...)) for v = Iterators.product(r...)]
            end
        end

              for k = keys(fields)])
    o = NamedTuple([k => (center + lb - o) / dx .+ 1 for (k, o) = pairs(field_origin)])
    _g = Dict([k => place(zeros(F, sizes[k]), g[k], o[k]) for k = keys(fields)])
    _center = round(center / dx) + 1 + common_left_pad_amount
    SourceInstance(f, keys(fields), g, _g, o, _center, label)
end

# Complex
apply(v::AbstractVector{<:SourceInstance}, t::Real; kw...) = apply(v, t, kw)
function apply(v::AbstractVector{<:SourceInstance}, t::Real, kw)
    dict([
        # sum([real(s.f(t) .* s._g[k]) for s = v if k in s.k], init=kw[k])
        k => begin

            a = [real(s.f(t) .* s._g[k]) for s = v if k in s.k]
            if isempty(a)
                kw[k]
            else
                kw[k] .+ sum(a)
            end
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