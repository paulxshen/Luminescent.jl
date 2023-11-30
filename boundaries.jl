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
    # function PML(dims, dx)
    #     new(dims, round(Int, 1 / dx))
    # end
end

struct Padding
    k
    b
    l
    r

end
function apply(b::AbstractVector{<:Padding}, fields,)
    for b = b
        fields = apply(b, fields,)
    end
    fields
end
function apply(p::Padding, fields::T,) where {T}
    @unpack l, r, b, k = p
    fields = merge(NamedTuple(fields), (; k => pad(PaddedArray, fields[k], b, l, r)))
    # ComponentArray(fields; (; k => pad(PaddedArray, fields[k], b, l, r))...)
    # ComponentArray(fields)
end
