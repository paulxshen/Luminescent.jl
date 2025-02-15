function getdimsperm(frame)
    r = zeros(Int, 3)
    vs = vcat(collect(eachcol(Matrix(I, 3, 3))), collect(eachcol(-Matrix(I, 3, 3))))
    for i = 1:3
        j = findfirst(vs) do v
            all(isapprox.(frame[:, i], v; atol=0.01))
        end
        isnothing(j) && return nothing
        if j > 3
            j = 3 - j
        end
        r[i] = j
    end
    r
end

function permutexyz(d, p, N=length(p))
    _p = @ignore_derivatives invperm(p)
    namedtuple([
        begin
            v = d[k]
            k = string(k)
            i = @ignore_derivatives findfirst(k[end], "xyz",)
            # global _a = p, i, k, d, _p
            k = @ignore_derivatives Symbol(k[1:end-1] * "xyz"[abs(_p[i])])
            k => sign(p[i]) * permutedims(v, p, N)
        end for k = keys(d)
    ])
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

function add_current_keys!(d::AbstractDict)
    for k in keys(d)
        if startswith(String(k), "E")
            s = String(k)[end]
            d[Symbol("J$s")] = deepcopy(d[k])
            d[Symbol("P$s")] = deepcopy(d[k])

        end
    end
    d
end

function groupkeys(d)
    r = dict()
    for k in keys(d)
        pre = String(k)[1] |> Symbol
        if !haskey(r, pre)
            r[pre] = dict()
        end
        r[pre][k] = d[k]
    end
    r
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
function cluster(v)
    groups = []
    group = []
    for (i, x) = enumerate(v)
        if isempty(group) || x / v[group[1]] < 1.1
            push!(group, i)
        else
            push!(groups, group)
            group = [i]
        end
    end
    if !isempty(group)
        push!(groups, group)
    end
    groups
end