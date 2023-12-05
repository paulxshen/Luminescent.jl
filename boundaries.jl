using LinearAlgebra, UnPack

struct Periodic
    dims
end

struct PEC
    dims
end
struct PMC
    dims
end
struct PEMC
    dims
end
struct PML
    dims
    d::Real
    σ::Real
    function PML(dims, d=0.5f0, σ=2.0f0)
        new(dims, d, σ)
    end
end

struct Padding
    k
    b
    l
    r
    info
    lazy
end
function apply(p::AbstractVector{<:Padding}, u, ;)
    u = NamedTuple(deepcopy(u))
    for p = p
        @unpack l, r, b, k, lazy, info = p
        if haskey(u, k)
            u = merge(u, (; k => begin
                a = u[k]
                y = pad(a, b, l, r; lazy)
                if info
                    y = PaddedArray(y, Int.(l .+ left(a)), Int.(r .+ right(a)))
                end
                y
            end))
        end
    end
    u
end