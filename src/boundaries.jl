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
    function PML(dims, d=0.5f0, σ=8.0f0)
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

function apply(p::AbstractVector{<:Padding}, a)
    y = a
    for p = p
        @unpack l, r, b = p
        y = pad(y, b, l, r;)
    end
    y
end
# function apply(p::Padding, a)
#     @unpack l, r, b = p
#     y = pad(a, b, l, r;)
# end
