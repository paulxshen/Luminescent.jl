using LinearAlgebra, UnPack
dl = Dict([
    -1 => "-x"
    1 => "+x"
    -2 => "-y"
    2 => "+y"
    3 => "+z"
    -3 => "-z"
])
"""
    Periodic(dims)

periodic boundary

"""
struct Periodic
    dims
end
Base.string(m::Periodic) = "Periodic on $(dl[m.dims]) side"
"""
PEC(dims)

perfect electrical conductor
dims: eg -1 for -x side
"""
struct PEC
    dims
end
Base.string(m::PEC) = "PEC ON $(dl[m.dims]) side"
"""
PMC(dims)

perfect magnetic conductor
"""
struct PMC
    dims
end
Base.string(m::PMC) = "PMC on $(dl[m.dims]) side"
struct PEMC
    dims
end
Base.string(m::PEMC) = "PEMC on $(dl[m.dims]) side"
"""
function PML(dims, d=0.25f0, σ=20.0f0)
    
    Constructs perfectly matched layers (PML aka ABC, RBC) boundary of depth `d` wavelengths 
    Doesn't need to be explictly declared as all unspecified boundaries default to PML
    """
struct PML
    dims
    d::Real
    σ::Real
    function PML(dims, d=0.5f0, σ=10.0f0)
        new(dims, d, σ)
    end
end
Base.string(m::PML) = "PML of depth $(m.d) on $(dl[m.dims]) side"

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
    default_size
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
    if length(p) == 1 && !isnothing(p[1].m)
        return a .* p[1].m
    end

    a_ = Buffer(a)
    # a_ .= a
    a_[axes(a)...] = a
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

function apply(v::AbstractVector{<:OutPad}, x::Real)
    # for p = v
    #     p.b != a && p.b != :replicate && error("cannot pad a constant with different values. it needs to be an array")
    # end
    # a
    # apply(p,)
    apply(v, x * ones(typeof(x), first(v).default_size))
end
function apply(p::AbstractVector{<:OutPad}, a)
    l = sum(getproperty.(p, :l))
    r = sum(getproperty.(p, :r))
    y = Buffer(a, Tuple(l .+ r .+ size(a)))
    y = place!(y, l .+ 1, a)
    for p = p
        l -= p.l
        r -= p.r
        y = pad!(y, p.b, p.l, p.r, l, r)
    end
    copy(y)
end
