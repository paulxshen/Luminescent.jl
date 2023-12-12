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


struct PlaneWave
    f
    C
    dims
end

struct GaussianBeam
    f
    σ
    C
    center
    dims
    function GaussianBeam(f, σ, center, dims; C...)
        new(f, σ, C, center, dims)
    end
end
struct CurrentSource
    f
    C
    lengths
    center
    function CurrentSource(f, lengths, center; C...)
        new(f, C, lengths, center)
    end
end

struct SourceEffect
    f
    g
    C
    start
end

function SourceEffect(s::PlaneWave, dx, sz, start)
    @unpack f, C, dims = s
    d = length(sz)
    g = ones([i == abs(dims) ? 1 : sz[i] for i = 1:d]...)
    start = start .+ dims < 0 ? 0 : [i == abs(dims) ? sz[i] - 1 : 0 for i = 1:d]
    SourceEffect(f, g, C, start)
end

function SourceEffect(s::GaussianBeam, dx, sz, start)
    @unpack f, σ, C, center, dims = s
    n = round(Int, 2σ / dx)
    r = n * dx
    I = [i == abs(dims) ? (0:0) : range(-r, r, length=(2n + 1)) for i = 1:length(center)]
    g = [gaussian(norm(F.(collect(v)))) for v = Iterators.product(I...)]
    start = start .- 1 .+ index(center, dx) .- round.(Int, (size(g) .- 1) ./ 2)
    SourceEffect(f, g, C, start)
end
function SourceEffect(s::CurrentSource, dx, sz, start)
    @unpack f, C, center, lengths = s
    n = round.(Int, lengths ./ dx)
    g = ones(n...)
    start = start .+ round.(Int, center ./ dx .- (n .- 1) ./ 2)
    SourceEffect(f, g, C, start)
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
                @unpack g, C, f, start = s
                if k in keys(C)
                    r = place(r, real(C[k] * f(t) .* g), start)
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
#         @unpack g, C, f, start = s
#         for f = keys(C)
#             u = merge(u, (; f => place(u[f], real(C[f] * s.f(t) .* g), start .+ left(u[f]))))
#         end
#     end
#     u
# end
