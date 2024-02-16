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

function SourceInstance(s::PlaneWave, dx, sizes, o, sz0)
    @unpack f, fields, dims, label = s
    d = length(first(sizes))
    g = Dict([k => fields[k] * ones([i == abs(dims) ? 1 : sz0[i] for i = 1:d]...) / dx for k = keys(fields)])
    o = NamedTuple([k =>
        o[k] .+ (dims < 0 ? 0 : [i == abs(dims) ? sizes[k][i] - 1 : 0 for i = 1:d])
                    for k = keys(o)])
    _g = Dict([k => place(zeros(F, sizes[k]), g[k], o[k]) for k = keys(fields)])
    c = NamedTuple([k => round.(Int, o[k] .+ size(first(values(g))) ./ 2) for k = keys(o)])
    SourceInstance(f, keys(fields), g, _g, o, c, label)
end

function SourceInstance(s::GaussianBeam, dx, sizes, o, stop)
    @unpack f, σ, fields, c, dims = s
    n = round(Int, 2σ / dx)
    r = n * dx
    I = [i == abs(dims) ? (0:0) : range(-r, r, length=(2n + 1)) for i = 1:length(c)]
    g = [gaussian(norm(F.(collect(v)))) for v = Iterators.product(I...)] / dx
    o = o .- 1 .+ index(c, dx) .- round.(Int, (size(g) .- 1) ./ 2)
    _g = place(zeros(F, sz), g, o)
    SourceInstance(f, keys(fields), g, _g, o, c, label)
end

function SourceInstance(s::Source, dx, sizes, o, stop)
    @unpack f, fields, c, lb, ub, label = s
    I = [a:dx:b for (a, b) = zip(lb, ub)]
    g = Dict([k => [fields[k](v...) for v = Iterators.product(I...)] / dx^count(lb .== ub) for k = keys(fields)])
    c = -1 .+ index(c, dx) .- round.(Int, (length.(I) .- 1) ./ 2)
    o = NamedTuple([k => o[k] .+ c for k = keys(o)])
    _g = Dict([k => place(zeros(F, sizes[k]), g[k], o[k]) for k = keys(fields)])
    c = NamedTuple([k => round.(Int, o[k] .+ size(first(values(g))) ./ 2) for k = keys(o)])
    SourceInstance(f, keys(fields), g, _g, o, c, label)
end

# function apply(s::SourceInstance, a, t::Real)
#     @unpack g, _g, f, o = s
#     # r = place(r, real(fields[k] * f(t) .* g), o)
#     a .+ real(f(t) .* _g)
# end
function apply(v::AbstractVector{<:SourceInstance}, t::Real; kw...)
    [
        sum([real(s.f(t) .* s._g[k]) for s = v if k in s.k], init=kw[k]) for k = keys(kw)
        # end for (k, a) = pairs(kw)
    ]
end

# function apply(d,t; kw...)
#     [apply(d[k], kw[k]) for k = keys(kw)]
# end