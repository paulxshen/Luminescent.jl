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
function apply(p::AbstractVector{<:Padding}, fields,)
    res = deepcopy(fields)
    for p = p
        @unpack l, r, b, k = p
        res = merge(NamedTuple(res), (; k => pad(PaddedArray, res[k], b, l, r)))
    end
    res
end