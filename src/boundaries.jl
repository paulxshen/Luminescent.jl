using LinearAlgebra, UnPack

"""
    Periodic(dims)

periodic boundary
"""
struct Periodic
    dims
end
"""
    PEC(dims)

perfect electrical conductor
"""
struct PEC
    dims
end
"""
    PMC(dims)

perfect magnetic conductor
"""
struct PMC
    dims
end
struct PEMC
    dims
end
"""
    function PML(dims, d=0.5f0, σ=8.0f0)

Constructs perfectly matched layers (PML aka ABC, RBC) boundary of depth `d` wavelengths 
"""
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

function apply(p; kw...)
    [apply(p[k], v) for (k, v) = pairs(kw)]
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
