getdimsperm(dims::Int) =
    if dims == 1
        [2, 3, 1]
    elseif dims == 2
        [-1, 3, 2]
    elseif dims == 3
        [1, 2, 3]
    end

function getdimsperm(L::Base.AbstractVecOrTuple)
    a = Int[]
    b = Int[]
    for (i, v) = enumerate(L)
        if abs(v) > 1e-3
            push!(a, sign(v) * i)
        else
            push!(b, i)
        end
    end
    v = vcat(a, b)
    v[end] *= sign(Permutation(abs.(v))) * prod(sign.(v))
    v
end
function permutexyz(d, p, N=length(p))
    _p = @ignore_derivatives invperm(p)
    namedtuple([
        begin
            k = string(k)
            i = findfirst(k[end], "xyz",)
            # global _a = p, i
            Symbol(k[1:end-1] * "xyz"[abs(p[i])]) => sign(p[i]) * permutedims(v, _p, N)
        end for (k, v) = pairs(d)
    ])
end

struct Grid
    F
    N
    L
    sz
    deltas
    lb
    field_lims
    field_sizes
    field_boundvals
    field_deltas
    field_diffdeltas
    field_diffpadvals
end
function get_polarization(u)
    # u=flatten(u)
    if haskey(u, :E)
        if length(u.E) == 2
            return :TE
        elseif length(u.H) == 2
            return :TM
        end
    else
        if haskey(u, :Hy)
            return :TE
        elseif haskey(u, :Ey)
            return :TM
        end
    end
    nothing
end
function add_current_keys(N)
    N = OrderedDict(pairs(N))
    add_current_keys!(N)
end
function add_current_keys!(N::AbstractDict)
    for k in keys(N) |> collect
        if startswith(String(k), "E")
            N[Symbol("J" * String(k)[2:end])] = N[k]
        end

        # if startswith(String(k), "H")
        #     N[Symbol("M" * String(k)[2:end])] = N[k]
        # end
    end
    N
end


_make_field_deltas(d::Real, a...) = d
function _make_field_deltas(d, N, field_boundvals, field_sizes, i, isdiff=false)
    NamedTuple([k => begin
        if isdiff * isnothing(v[i, 1])
            d = vcat(d[1], (d[1:end-1] .+ d[2:end]) ./ 2)
        end
        sel = i .== 1:N
        reshape(d, Tuple(1 - sel + field_sizes[k][i] * sel))
    end for (k, v) = pairs(field_boundvals)])
end
