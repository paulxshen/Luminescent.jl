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
dims: eg -1 for -x side
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
    function PML(dims, d=0.25f0, σ=20.0f0)

Constructs perfectly matched layers (PML aka ABC, RBC) boundary of depth `d` wavelengths 
Doesn't need to be explictly declared as all unspecified boundaries default to PML
"""
struct PML
    dims
    d::Real
    σ::Real
    function PML(dims, d=0.25f0, σ=20.0f0)
        new(dims, d, σ)
    end
end

struct Padding

    b
    l
    r
    out
    # info
    # lazy
end


function apply(p::AbstractVector{<:Padding}, a)
    y = a
    for p = p
        @unpack l, r, b, out = p
        y = out ? pad(y, b, l, r;) : pad!(y, b, l, r;)
    end
    y
end
# function apply(p::Padding, a)
#     @unpack l, r, b = p
#     y = pad(a, b, l, r;)
# end
