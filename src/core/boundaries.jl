using LinearAlgebra, UnPack
dl = OrderedDict([
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
    m
    ramp_frac
    function PML(dims, d, σ, m; ramp_frac=0.2f0)
        new(dims, d, σ, m, ramp_frac)
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
left(x::InPad) = x.l
right(x::InPad) = x.r

struct OutPad
    b
    l
    r
    default_size
end
# @functor OutPad


# function apply!(p::AbstractVector{<:InPad}, a::AbstractArray; nonzero_only=false)
#     if nonzero_only
#         p = filter(p) do p
#             p.b != 0
#         end
#     end
#     isempty(p) && return a

#     for p = p
#         @unpack l, r, b, m = p
#         if isnothing(m)
#             pad!(a, b, l, r;)
#         else
#             a .= a .* m
#         end
#     end
#     a
# end
