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
    function PlaneWave(f, dims; mode...)

Constructs plane wave source

Args
- f: time function
- dims: eg -1 for wave coming from -x face
- mode: which mode to excite & their scaling constants (typically a current source, eg Jz=1)
"""
struct PlaneWave
    f
    mode
    dims
    label
    function PlaneWave(f, dims, label=""; mode...)
        new(f, mode, dims, label)
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
    meta
    function Source(f, mode, c, lb, ub, label::String=""; kw...)
        new(f, mode, c, lb, ub, label, kw)
    end
    function Source(f, mode, c, L, label::String=""; kw...)
        new(f, mode, c, -L / 2, L / 2, label, kw)
    end
end
# mode(m::Source) = m.mode
Base.string(m::Union{Source,ModalSource}) =
    """
    $(m.label): $(count((m.ub.-m.lb).!=0))-dimensional source, centered at $(m.center|>d2), spanning from $(m.lb|>d2) to $(m.ub|>d2) relative to center,"""
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
    meta
end
@functor SourceInstance (g, _g,)
Porcupine.keys(m::SourceInstance) = keys(m.g)

function SourceInstance(s::ModalSource, Δ, sizes, origin, fieldlims, common_left_pad_amount, ; F=Float32)
    @unpack f, center, lb, ub, normal, tangent, meta = s
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
    sz = max.(1, round.(Int, L ./ Δ[1:end-1])) |> Tuple
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
        zaxis = convert.(F, normal)
        xaxis = convert.(F, tangent)
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
        mode = (; Jx, Jy)
    else
        mode = (; Jx, Jy, Jz)
    end
    v = zip(lb, ub)
    lb = minimum.(v)
    ub = maximum.(v)
    mode = ignore_derivatives() do
        mode / Δ[findfirst(abs.(zaxis) .> 0.001)]
    end
    SourceInstance(Source(f, mode, center, lb, ub; meta...), Δ, sizes, origin, fieldlims, common_left_pad_amount; F)
end
function SourceInstance(s::PlaneWave, Δ, sizes, common_left_pad_amount, fl, sz0; F=Float32)
    @unpack f, mode, dims, label, meta = s
    _F(x::Real) = F(x)
    _F(x::Complex) = ComplexF32(x)
    f = _F ∘ f
    d = length(common_left_pad_amount)
    g = Dict([k => _F.(mode[k]) * ones([i == abs(dims) ? 1 : sz0[i] for i = 1:d]...) / Δ for k = keys(mode)])
    o = NamedTuple([k =>
        1 .+ fl[k] .+ (dims < 0 ? 0 : [i == abs(dims) ? sizes[k][i] - 1 : 0 for i = 1:d])
                    for k = keys(fl)])
    _g = Dict([k => place(zeros(F, sizes[k]), o[k], g[k],) for k = keys(mode)])
    c = first(values(sizes)) .÷ 2
    SourceInstance(f, g, _g, o, c, label, meta)
end

function SourceInstance(s::GaussianBeam, Δ, sizes, fl, stop; F=Float32)
    _F(x::Real) = F(x)
    _F(x::Complex) = ComplexF32(x)
    f = _F ∘ f
    @unpack f, σ, mode, c, dims = s
    n = round(Int, 2σ / Δ)
    r = n * Δ
    r = [i == abs(dims) ? (0:0) : range(-r, r, length=(2n + 1)) for i = 1:length(c)]
    g = [gaussian(norm(F.(collect(v)))) for v = Iterators.product(r...)] / Δ
    fl = fl .- 1 .+ index(c, Δ) .- round.(Int, (size(g) .- 1) ./ 2)
    _g = place(zeros(F, sz), g, fl)
    SourceInstance(f, g, _g, fl, c, label)
end

function SourceInstance(s::Source, Δ, sizes, origin, fieldlims, common_left_pad_amount; F=Float32)
    @unpack f, mode, center, lb, ub, label, meta = s
    _F(x::Real) = F(x)
    _F(x::Complex) = ComplexF32(x)

    f = _F ∘ f
    # Δmode = Δ[end-N+2:end]
    # g = Dict([k =>
    #     begin
    #         if isa(mode[k], AbstractArray)
    #             sz0 = size(mode[k])
    #             sz = max.(1, round.(Int, abs.(ub - lb) ./ Δ)) |> Tuple
    #             # sz0 != sz && @warn "source array size$sz0 not same  as domain size $sz. source will be interpolated"
    #             imresize(_F.(mode[k]), sz, method=ImageTransformations.Lanczos4OpenCV())
    #         else
    #             # r = [a == b ? (a:a) : (a+Δ/2*sign(b - a):Δ*sign(b - a):b) for (a, b) = zip(lb, ub)]
    #             # [_F.(mode[k](v...)) for v = Iterators.product(r...)]
    #         end
    #     end for k = keys(mode)])
    g = mode
    o = NamedTuple([k => F.((center + lb - origin) ./ Δ - fl[:, 1] .+ 1.5) for (k, fl) = pairs(fieldlims)])
    @show center, lb, origin, Δ, fieldlims, o
    global _g = Dict([k => begin
        a = zeros(ComplexF32, sizes[k])
        setindexf!(a, g[k], range.(o[k], o[k] + size(g[k]) - 1)...)
        a
    end for k = keys(mode)])
    # error("stop")
    _center = round.(Int, center ./ Δ) + 1 + common_left_pad_amount
    SourceInstance(f, g, _g, o, _center, label, meta)
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