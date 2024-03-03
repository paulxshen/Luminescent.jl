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
    c
    lb
    ub
    label
    function Source(f, c, lb, ub, label::String=""; fields...)
        new(f, fields, c, lb, ub, label)
    end
end

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
@functor SourceInstance

function SourceInstance(s::PlaneWave, dx, sizes, lc, fl, sz0; F=Float32)
    @unpack f, fields, dims, label = s
    f = _F ∘ f
    d = length(lc)
    g = Dict([k => _F(fields[k]) * ones([i == abs(dims) ? 1 : sz0[i] for i = 1:d]...) / dx for k = keys(fields)])
    o = NamedTuple([k =>
        1 .+ fl[k] .+ (dims < 0 ? 0 : [i == abs(dims) ? sizes[k][i] - 1 : 0 for i = 1:d])
                    for k = keys(fl)])
    _g = Dict([k => place(zeros(F, sizes[k]), g[k], o[k]) for k = keys(fields)])
    c = first(values(sizes)) .÷ 2
    SourceInstance(f, keys(fields), g, _g, o, c, label)
end

function SourceInstance(s::GaussianBeam, dx, sizes, fl, stop; F=Float32)
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

function SourceInstance(s::Source, dx, sizes, lc, fl, stop; F=Float32)
    @unpack f, fields, c, lb, ub, label = s
    # println(fl)

    f = _F ∘ f
    r = [vcat(a:dx:0, dx:dx:b) for (a, b) = zip(lb, ub)]
    C = 1 /
        dx^count(lb .== ub)
    g = Dict([k =>
        C * begin
            if isa(fields[k], AbstractArray)
                # imresize(fields[k], ratio=1)
                imresize(_F.(fields[k]), Tuple(length.(r)), method=ImageTransformations.Lanczos4OpenCV())
            else
                [_F.(fields[k](v...)) for v = Iterators.product(r...)]
            end
        end

              for k = keys(fields)])
    c = round.(c ./ dx) .+ 1
    o = c + round.(lb ./ dx)
    c = c + lc
    o = NamedTuple([k => fl[k] .+ o for k = keys(fl)])
    _g = Dict([k => place(zeros(F, sizes[k]), g[k], o[k]) for k = keys(fields)])
    SourceInstance(f, keys(fields), g, _g, o, c, label)
end

_F(x::Real) = F(x)
_F(x::Complex) = complex(F)(x)
# Complex
function apply(v::AbstractVector{<:SourceInstance}, t::Real; kw...)
    [
        # sum([real(s.f(t) .* s._g[k]) for s = v if k in s.k], init=kw[k])
        begin

            a = [real(s.f(t) .* s._g[k]) for s = v if k in s.k]
            if isempty(a)
                kw[k]
            else
                kw[k] .+ sum(a)
            end
        end
        for k = keys(kw)
        # end for (k, a) = pairs(kw)
    ]
end

# function apply(d,t; kw...)
#     [apply(d[k], kw[k]) for k = keys(kw)]
# end