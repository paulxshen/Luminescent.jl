using LinearAlgebra, UnPack, NamedTupleTools
# ArrayPadding,
include("../ArrayPadding.jl/src/pad.jl")
# __precompile__(false)

gaussian(x; μ=0, σ=1) = exp(-((x - μ) / σ)^2)

function place(a, b, start; lazy=false)
    # @show size(a), size(b), start
    a + pad(b, 0, Tuple(start) .- 1, size(a) .- size(b) .- Tuple(start) .+ 1; lazy)
end

function place(a, b; center)
    # @show size(a), size(b), center
    place(a, b, center .- floor.((size(b) .- 1) .÷ 2))
end

function fnt(a::NamedTuple)
    r = (;)
    ignore() do
        for (k, v) = pairs(a)
            f, c = String(k)
            f, c = Symbol.([f, c])
            r = merge_recursive(r, (; f => (; c => v)))
        end
    end
    r
    # NamedTuple([s[1] => (; s[2] => 1 for )) for s = String.(a)])

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
    start = start .+ round.(Int, center ./ dx) .- n
    SourceEffect(f, g, C, start)
end
function SourceEffect(s::CurrentSource, dx, sz, start)
    @unpack f, C, center, lengths = s
    n = round.(Int, lengths ./ dx)
    g = ones(n...)
    start = start .+ round.(Int, center ./ dx .- (n .- 1) ./ 2)
    SourceEffect(f, g, C, start)
end

function apply(s::AbstractVector{<:SourceEffect}, u, t)
    u = NamedTuple(deepcopy(u))
    for s = s
        @unpack g, C, f, start = s
        for f = keys(C)
            u = merge(u, (; f => place(u[f], real(C[f] * s.f(t) .* g), start .+ left(u[f]))))
        end
    end
    u
end
