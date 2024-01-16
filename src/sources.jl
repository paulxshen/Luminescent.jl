# __precompile__(false)
(m::Number)(a...) = m
gaussian(x; μ=0, σ=1) = exp(-((x - μ) / σ)^2)

function place(a, b, start; lazy=false)
    a + pad(b, 0, Tuple(start) .- 1, size(a) .- size(b) .- Tuple(start) .+ 1; lazy)
    # buf = Buffer(a)
    # for i = eachindex(a)
    #     buf[i] = a[i]
    # end
    # buf[[i:j for (i, j) = zip(start, start .+ size(b) .- 1)]...] += b
    # copy(buf)
    # stop = start .+ size(b) .- 1
    # map(Iterators.product(axes(a)...)) do I
    #     all(start .<= I .<= stop) ? b[(I .- start .+ 1)...] + a[I...] : a[I...]
    # end
end

function place(a, b; center)
    # @show size(a), size(b), center
    place(a, b, center .- floor.((size(b) .- 1) .÷ 2))
end



"""
    function PlaneWave(f, dims; fields...)

Constructs plane wave source

Args
- f: time function
- dims: eg -1 for wave coming from -x edge
- fields: which fields to excite & their scaling constants (typically a current source, eg Jz=1)
"""
struct PlaneWave
    f
    fields
    dims
    function PlaneWave(f, dims; fields...)
        new(f, fields, dims)
    end
end

"""
    function GaussianBeam(f, σ, center, dims; fields...)

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
    center
    dims
    function GaussianBeam(f, σ, center, dims; fields...)
        new(f, σ, fields, center, dims)
    end
end
"""
    function Source(f, center, bounds; fields...)
    function Source(f, center, L::AbstractVector{<:Real}; fields...)

Constructs custom centered source. Can be used to specify modal sources

Args
- f: time function
- L: source dimensions in [wavelengths]
- fields: which fields to excite & their scaling constants (typically a current source, eg Jz=1)
"""
struct Source
    f
    fields
    center
    bounds
    function Source(f, center, bounds; fields...)
        new(f, fields, center, bounds)
    end
end
function Source(f, center, L::Union{AbstractVector{<:Real},Tuple{<:Real}}; fields...)
    Source(f, center, [[a, a] for a = L]; fields...)
end
"""
    function UniformSource(f, center, L; fields...)

Constructs uniform source

Args
- f: time function
- L::Vector: lengths in [wavelengths]
- fields: which fields to excite & their scaling constants (typically a current source, eg Jz=1)
"""
struct UniformSource
    f
    fields
    L
    center
    function UniformSource(f, center, L; fields...)
        new(f, fields, L, center)
    end
end

struct SourceEffect
    f
    g
    _g
    fields
    start
end

function SourceEffect(s::PlaneWave, dx, sz, start, stop)
    @unpack f, fields, dims = s
    d = length(sz)
    a = ones([i == abs(dims) ? 1 : sz[i] for i = 1:d]...) / dx
    start = start .+ (dims < 0 ? 0 : [i == abs(dims) ? sz[i] - 1 : 0 for i = 1:d])
    g = Dict([k => fields[k] * a for k = keys(fields)])
    _g = Dict([k => place(zeros(F, sz), g[k], start) for k = keys(fields)])
    SourceEffect(f, g, _g, fields, start)
end

function SourceEffect(s::GaussianBeam, dx, sz, start, stop)
    @unpack f, σ, fields, center, dims = s
    n = round(Int, 2σ / dx)
    r = n * dx
    I = [i == abs(dims) ? (0:0) : range(-r, r, length=(2n + 1)) for i = 1:length(center)]
    g = [gaussian(norm(F.(collect(v)))) for v = Iterators.product(I...)] / dx
    start = start .- 1 .+ index(center, dx) .- round.(Int, (size(g) .- 1) ./ 2)
    _g = place(zeros(F, sz), g, start)
    SourceEffect(f, g, _g, fields, start)
end
function SourceEffect(s::Source, dx, sz, start, stop)
    @unpack f, fields, center, bounds = s
    # R = round.(Int, L ./ 2 / dx)
    # I = range.(-R, R)
    # I = [round(Int, a / dx):round(Int, b / dx) for (a, b) = bounds]
    bounds = [isa(b, Number) ? [b, b] : b for b = bounds]
    I = [b[1]:dx:b[2] for b = bounds]
    g = Dict([k => [fields[k](v...) for v = Iterators.product(I...)] / dx^count(getindex.(bounds, 1) .== getindex.(bounds, 2)) for k = keys(fields)])
    start = start .- 1 .+ index(center, dx) .- round.(Int, (length.(I) .- 1) ./ 2)
    _g = Dict([k => place(zeros(F, sz), g[k], start) for k = keys(fields)])
    SourceEffect(f, g, _g, fields, start)
    # n = max.(1, round.(Int, L ./ dx))
    # g = ones(n...) / dx^count(L .== 0)
    # start = start .+ round.(Int, center ./ dx .- (n .- 1) ./ 2)
end

function apply(s::AbstractVector{<:SourceEffect}, t; kw...)
    k = 0
    ignore() do
        k = keys(kw)
    end
    [
        begin
            r = kw[k]
            for s = s
                @unpack g, _g, fields, f, start = s
                if k in keys(fields)
                    # r = place(r, real(fields[k] * f(t) .* g), start)
                    r = r .+ real(f(t) .* _g[k])
                end
            end
            r
        end for k = k
        # end for (k, a) = pairs(kw)
    ]
end
# function apply(s::AbstractVector{<:SourceEffect}, u, t)
#     u = NamedTuple(deepcopy(u))
#     for s = s
#         @unpack g, fields, f, start = s
#         for f = keys(fields)
#             u = merge(u, (; f => place(u[f], real(fields[f] * s.f(t) .* g), start .+ left(u[f]))))
#         end
#     end
#     u
# end
