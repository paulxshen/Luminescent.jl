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
    function PML(dims, d=0.25f0, σ=10.0f0)
        new(dims, d, σ)
    end
end

struct InPad

    b
    l
    r
    m
    function InPad(b, l, r, m=nothing)
        new(b, l, r, m)
    end
end
@functor InPad (m,)
struct OutPad
    b
    l
    r
end
# @functor OutPad


function apply!(p::AbstractVector{<:InPad}, a::AbstractArray)
    for p = p
        @unpack l, r, b, m = p
        if isnothing(m)
            pad!(a, b, l, r;)
        else
            a .= a .* m
        end
    end
    a
end
function apply(p::AbstractVector{<:InPad}, a::AbstractArray)
    a_ = Buffer(a)
    for p = p
        @unpack l, r, b, m = p
        if isnothing(m)
            pad!(a_, b, l, r;)
        else
            a_[axes(a)...] = a .* m
        end
    end

    copy(a_)
end

function apply(p::AbstractVector{<:OutPad}, a)
    l = sum(getproperty.(p, :l))
    r = sum(getproperty.(p, :r))
    y = Buffer(a, Tuple(l .+ r .+ size(a)))
    y = place!(y, a, l .+ 1)
    for p = p
        l -= p.l
        r -= p.r
        y = pad!(y, p.b, p.l, p.r, l, r)
    end
    copy(y)
end
# function apply(p::InPad, a)
#     @unpack l, r, b = p
#     y = pad(a, b, l, r;)
# end
