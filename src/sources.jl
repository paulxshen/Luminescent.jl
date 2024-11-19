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
    tags
    function ModalSource(f, mode, center::Base.AbstractVecOrTuple, normal, tangent, lb, ub, ; tags=Dict())
        new(f, center, lb, ub, normal, tangent, E2J(mode), string(label), tags)
    end
    function ModalSource(f, mode, center::Base.AbstractVecOrTuple, normal, tangent, L, ; tags=Dict())
        new(f, center, -L / 2, L / 2, normal, tangent, E2J(mode), string(label), tags)
    end

end
"""
    function PlaneWave(f, dims; mode...)

Constructs plane wave source

Args
- f: time function
- dims: eg -1 for wave coming from -x face
- mode: which mode to excite & their scaling constants (typically a current source, eg Jz=1)
"""
struct PlaneWave
    λmodes
    dims
    label
    function PlaneWave(λmodes, dims, label="")
        new(λmodes, dims, label)
    end
end
@functor PlaneWave
"""
    function GaussianBeam(f, σ, c, dims; mode...)

Constructs gaussian beam source

Args
- f: time function
- σ: std dev length
- dims: eg 1 for x direction
- mode: which mode to excite & their scaling constants (typically a current source, eg Jz=1)
"""
struct GaussianBeam
    f
    σ
    mode
    c
    dims
    function GaussianBeam(f, σ, c, dims; mode...)
        new(f, σ, mode, c, dims)
    end
end
@functor GaussianBeam

"""
    function Source(f, c, lb, ub, label=""; mode...)
    function Source(f, c, L, label=""; mode...)
        
Constructs custom  source. Can be used to specify uniform or modal sources

Args
- f: time function
- c: origin or center of source
- lb: lower bounds wrt to c
- ub: upper bounds wrt to c
- L: source dimensions in [wavelengths]
- mode: which mode to excite & their scaling constants (typically a current source, eg Jz=1)
"""
struct Source
    f
    mode
    center
    lb
    ub
    label
    tags
    function Source(f, mode, c, lb, ub; kw...)
        new(f, mode, c, lb, ub, string(label), kw)
    end
    function Source(f, mode, c, L; kw...)
        new(f, mode, c, -L / 2, L / 2, string(label), kw)
    end
end
# mode(m::Source) = m.mode
# Base.string(m::Union{Source,ModalSource}) =
#     """
#     $(m.label): $(count((m.ub.-m.lb).!=0))-dimensional source, centered at $(m.center|>d2), spanning from $(m.lb|>d2) to $(m.ub|>d2) relative to center,"""
#  exciting $(join(keys(m.mode),", "))"""

function Source(f, c, L, label::AbstractString=""; mode...)
    # Source
    # function Source(f, c, L::Union{AbstractVector{<:Real},Tuple{<:Real}}; mode...)
    Source(f, c, -L / 2, L / 2, label; mode...)
end
Source = Source

struct SourceInstance
    f
    g
    _g
    o
    center
    label
    tags
end
@functor SourceInstance (g, _g,)
Porcupine.keys(m::SourceInstance) = keys(m.g)

function SourceInstance(s::ModalSource, deltas, sizes, origin, field_lims, ; F=Float32)
    @unpack f, center, lb, ub, normal, tangent, label, tags = s
    J = OrderedDict()
    for k = (:Jx, :Jy, :Jz)
        if k in keys(s.mode)
            J[k] = ComplexF32.(s.mode[k])
        else
            J[k] = zeros(ComplexF32, size(s.mode(1)))
        end
    end
    J = values(J)

    L = ub .- lb
    d = length(lb)
    D = length(center) # 2D or 3D
    if D == 2
        zaxis = [normal..., 0]
        yaxis = [0, 0, 1]
        # xaxis = cross(yaxis, zaxis)
        xaxis = cross(yaxis, zaxis)
    else
        zaxis = convert.(F, normal)
        xaxis = convert.(F, tangent)
        yaxis = cross(zaxis, xaxis)
    end

    frame = [xaxis, yaxis, zaxis]
    J = reframe(frame, J)
    if D == 2
        J = [J[:, :, 1] for J in J]
    end
    lb = sum(lb .* frame[1:d])[1:D]
    ub = sum(ub .* frame[1:d])[1:D]
    L = ub - lb

    Jx, Jy, Jz = J
    if D == 2
        mode = (; Jx, Jy)
    else
        mode = (; Jx, Jy, Jz)
    end
    v = zip(lb, ub)
    lb = minimum.(v)
    ub = maximum.(v)
    n = findfirst(abs.(zaxis) .> 0.001)
    mode = ignore_derivatives() do
        mode / deltas[n][1]#[findfirst(>(center[n]), cumsum(deltas[n]))]
    end
    SourceInstance(Source(f, mode, center, lb, ub; label, tags...), deltas, sizes, origin, field_lims, ; F)
end

function SourceInstance(s::PlaneWave, deltas, sizes, origin, field_lims, ; F=Float32)
    @unpack f, center, lb, ub, normal, tangent, label, tags = s
    J = OrderedDict()
    for k = (:Jx, :Jy, :Jz)
        if k in keys(s.mode)
            J[k] = ComplexF32.(s.mode[k])
        else
            J[k] = zeros(ComplexF32, size(s.mode(1)))
        end
    end
    J = values(J)

    L = ub .- lb
    d = length(lb)
    D = length(center) # 2D or 3D
    if D == 2
        zaxis = [normal..., 0]
        yaxis = [0, 0, 1]
        # xaxis = cross(yaxis, zaxis)
        xaxis = cross(yaxis, zaxis)
    else
        zaxis = convert.(F, normal)
        xaxis = convert.(F, tangent)
        yaxis = cross(zaxis, xaxis)
    end

    frame = [xaxis, yaxis, zaxis]
    J = reframe(frame, J)
    if D == 2
        J = [J[:, :, 1] for J in J]
    end
    lb = sum(lb .* frame[1:d])[1:D]
    ub = sum(ub .* frame[1:d])[1:D]
    L = ub - lb

    Jx, Jy, Jz = J
    if D == 2
        mode = (; Jx, Jy)
    else
        mode = (; Jx, Jy, Jz)
    end
    v = zip(lb, ub)
    lb = minimum.(v)
    ub = maximum.(v)
    n = findfirst(abs.(zaxis) .> 0.001)
    mode = ignore_derivatives() do
        mode / deltas[n][1]#[findfirst(>(center[n]), cumsum(deltas[n]))]
    end
    SourceInstance(Source(f, mode, center, lb, ub; label, tags...), deltas, sizes, origin, field_lims, ; F)
end

function SourceInstance(s::Source, deltas, sizes, origin, field_lims, ; F=Float32)
    @unpack f, mode, center, lb, ub, label, tags = s
    _F(x::Real) = F(x)
    _F(x::Complex) = ComplexF32(x)

    f = _F ∘ f
    start = v2i(center + lb - origin, deltas)
    stop = v2i(center + ub - origin, deltas)
    sel = abs.(stop - start) .>= 0.001
    start += 0.5sel
    stop -= 0.5sel

    sz0 = size(mode(1))
    sz = stop - start .+ 1
    @show sz0, sz
    # g = NamedTuple([k => imresize(a, Tuple(round(sz))) for (k, a) = pairs(mode)])
    g = mode

    o = NamedTuple([k => F.(start - fl[:, 1] + 1) for (k, fl) = pairs(field_lims)])
    # @show center, lb, origin, deltas, field_lims, o
    _g = namedtuple([k => begin
        a = zeros(ComplexF32, sizes[k])
        setindexf!(a, g[k], range.(o[k], o[k] + size(g[k]) - 1)...)
        a
    end for k = keys(mode)])
    # error("stop")
    _center = round(v2i(center - origin, deltas) + 0.5)
    @show center, origin, _center

    SourceInstance(f, g, _g, o, _center, label, tags)
end

# Complex
function inject_sources(v::AbstractVector{<:SourceInstance}, t::Real)
    ks = keys(v[1])
    namedtuple([
        k => sum([real(s.f(t) .* s._g[k]) for s = v if k in keys(s)])
        for k = ks])
end

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


# function SourceInstance(s::PlaneWave, deltas, sizes, common_left_pad_amount, fl, sz0; F=Float32)
#     @unpack f, mode, dims, label, tags = s
#     _F(x::Real) = F(x)
#     _F(x::Complex) = ComplexF32(x)
#     f = _F ∘ f
#     d = length(common_left_pad_amount)
#     g = Dict([k => _F.(mode[k]) * ones([i == abs(dims) ? 1 : sz0[i] for i = 1:d]...) / deltas for k = keys(mode)])
#     o = NamedTuple([k =>
#         1 .+ fl[k] .+ (dims < 0 ? 0 : [i == abs(dims) ? sizes[k][i] - 1 : 0 for i = 1:d])
#                     for k = keys(fl)])
#     _g = Dict([k => place(zeros(F, sizes[k]), o[k], g[k],) for k = keys(mode)])
#     c = first(values(sizes)) .÷ 2
#     SourceInstance(f, g, _g, o, c, label, tags)
# end

# function SourceInstance(s::GaussianBeam, deltas, sizes, fl, stop; F=Float32)
#     _F(x::Real) = F(x)
#     _F(x::Complex) = ComplexF32(x)
#     f = _F ∘ f
#     @unpack f, σ, mode, c, dims = s
#     n = round(Int, 2σ / deltas)
#     r = n * deltas
#     r = [i == abs(dims) ? (0:0) : range(-r, r, length=(2n + 1)) for i = 1:length(c)]
#     g = [gaussian(norm(F.(collect(v)))) for v = Iterators.product(r...)] / deltas
#     fl = fl .- 1 .+ index(c, deltas) .- round.(Int, (size(g) .- 1) ./ 2)
#     _g = place(zeros(F, sz), g, fl)
#     SourceInstance(f, g, _g, fl, c, label)
# end
