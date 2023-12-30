# __precompile__(false)

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
    function CenteredSource(f, g, center, L; fields...)

Constructs custom centered source. Can be used to specify modal sources

Args
- f: time function
- g: spatial function
- L: source dimensions in [wavelengths]
- fields: which fields to excite & their scaling constants (typically a current source, eg Jz=1)
"""
struct CenteredSource
    f
    g
    fields
    center
    L
    function CenteredSource(f, g, center, L; fields...)
        new(f, g, fields, center, L)
    end
end

"""
    function UniformSource(f, lengths, center; fields...)

Constructs uniform source

Args
- f: time function
- fields: which fields to excite & their scaling constants (typically a current source, eg Jz=1)
"""
struct UniformSource
    f
    fields
    lengths
    center
    function UniformSource(f, lengths, center; fields...)
        new(f, fields, lengths, center)
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
    g = ones([i == abs(dims) ? 1 : sz[i] for i = 1:d]...)
    start = start .+ (dims < 0 ? 0 : [i == abs(dims) ? sz[i] - 1 : 0 for i = 1:d])
    _g = place(zeros(F, sz), g, start)
    SourceEffect(f, g, _g, fields, start)
end

function SourceEffect(s::GaussianBeam, dx, sz, start, stop)
    @unpack f, σ, fields, center, dims = s
    n = round(Int, 2σ / dx)
    r = n * dx
    I = [i == abs(dims) ? (0:0) : range(-r, r, length=(2n + 1)) for i = 1:length(center)]
    g = [gaussian(norm(F.(collect(v)))) for v = Iterators.product(I...)]
    start = start .- 1 .+ index(center, dx) .- round.(Int, (size(g) .- 1) ./ 2)
    _g = place(zeros(F, sz), g, start)
    SourceEffect(f, g, _g, fields, start)
end
function SourceEffect(s::CenteredSource, dx, sz, start, stop)
    @unpack f, g, fields, center, L = s
    R = round.(Int, L ./ 2 / dx)
    I = range.(-R, R)
    g = [g(dx .* v...) for v = Iterators.product(I...)]
    start = start .- 1 .+ index(center, dx) .- round.(Int, (size(g) .- 1) ./ 2)
    _g = place(zeros(F, sz), g, start)
    SourceEffect(f, g, _g, fields, start)
end
function SourceEffect(s::UniformSource, dx, sz, start, stop)
    @unpack f, fields, center, lengths = s
    n = round.(Int, lengths ./ dx)
    g = ones(n...)
    start = start .+ round.(Int, center ./ dx .- (n .- 1) ./ 2)
    _g = place(zeros(F, sz), g, start)
    SourceEffect(f, g, _g, fields, start)
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
                    r = r .+ real(fields[k] * f(t) .* _g)
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
