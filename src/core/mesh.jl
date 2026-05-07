struct RangeFinder
    lbs
    ubs
    lbixs
    ubixs
    function RangeFinder(bboxes)
        v = sort.(collect.(enumerate.(eachrow(bboxes[:, 1, :]))), by=x -> x[2])
        lbs = reduce(hcat, [last.(v) for v = v])'
        lbixs = reduce(hcat, [first.(v) for v = v])'
        v = sort.(collect.(enumerate.(eachrow(bboxes[:, 2, :]))), by=x -> x[2])
        ubs = reduce(hcat, [last.(v) for v = v])'
        ubixs = reduce(hcat, [first.(v) for v = v])'
        new(lbs, ubs, lbixs, ubixs)
    end
end
function Base.findall(u, v, c::RangeFinder)
    @unpack lbs, ubs, lbixs, ubixs = c
    n = size(lbs, 2)
    r = 1:n
    for (u, v, lbs, ubs, lbixs, ubixs) = zip(u, v, eachrow(lbs), eachrow(ubs), eachrow(lbixs), eachrow(ubixs))
        i = searchsortedlast(lbs, v)
        i == 0 && return Int[]
        s = lbixs[1:i]
        r = intersect(r, s)
        isempty(r) && return Int[]

        j = searchsortedfirst(ubs, u)
        j == n + 1 && return Int[]
        s = ubixs[j:end]
        r = intersect(r, s)
        isempty(r) && return Int[]
    end
    r
end
Base.findall(v, c::RangeFinder) = Base.findall(v, v, c)

function share(v, x)
    v + x * v .^ 2 / sum(v .^ 2)
end

const ℛ = 1.4f0
function partition(l::T, v, dl, dr) where T<:Real
    n = ceil(Int, l * v)
    # @show l, v, n
    d0 = l / n
    if dl === dr === nothing
        return fill(d0, n)
    end
    if isnothing(dl)
        return reverse(partition(l, v, dr, dl))
    end
    if isnothing(dr)
        l < dl && return [l]
        ds = [dl]
        r = l - dl
        while true
            d = min(d0, ℛ * ds[end]) |> T
            if r > d
                push!(ds, d)
                r -= d
            elseif r > d / ℛ
                push!(ds, r)
                return ds
            else
                return share(ds, r)
            end
        end
    end

    xl = l * dl / (dl + dr)
    xr = l * dr / (dl + dr)
    (xl < dl / ℛ || xr < dr / ℛ) && return [l]
    vcat(partition(xl, v, dl, nothing), partition(xr, v, nothing, dr))
end

function makemesh(mvs, bbox, nres; bg, tol)
    @debug (; mvs, bbox)
    xs = Any[Any[a, b] for (a, b) = eachrow(bbox)]
    vs = Any[Any[bg] for _ = 1:length(xs)]

    for (m, v) = reverse(collect(mvs))
        x = boundingbox(m)
        a = ustrip.(getfield(coords(x.min), :coords))
        b = ustrip.(getfield(coords(x.max), :coords))
        for (dim, (xs, vs, a, b)) = enumerate(zip(xs, vs, a, b))
            a = max(xs[1], a)
            b = min(xs[end], b)
            if b > a
                i = searchsortedfirst(xs, a; by=first)
                j = searchsortedfirst(xs, b; by=first)
                insert!(xs, j, b)
                insert!(vs, j, vs[j-1](dim))
                if i < j
                    # vs[i:j-1] .= v
                    for k = i:j-1
                        vs[k] = max(vs[k](dim), v(dim))
                    end
                end
                insert!(xs, i, a)
                insert!(vs, i, v(dim))
            end
        end
    end

    for (xs, vs) = zip(xs, vs)
        c = 0
        for (i, d) = enumerate(diff(xs))
            if d < tol
                deleteat!(xs, i - c)
                deleteat!(vs, i - c)
                c += 1
            end
        end
    end

    vs *= nres
    _vs = deepcopy(vs)
    for (i, vs) = enumerate(vs)
        for (j, v) = enumerate(vs)
            if 1 < j < length(vs)
                vl = vs[j-1]
                vr = vs[j+1]
                if vl > v < vr
                    _vs[i][j] = 2max(vl, vr)
                end
            end
        end
    end

    ds = diff.(xs)
    Δs = [Any[nothing for _ = 1:length(ds)] for ds = ds]
    for (ds, _vs, Δs) = zip(ds, _vs, Δs)
        n = length(ds)
        for i = sort(1:n, by=i -> ds[i] / ceil(ds[i] * _vs[i]))
            # if i == 1
            #     dl = vl = nothing
            # else
            #     dl = Δs[i-1][end]
            # end
            # if i == n
            #     dr = vr = nothing
            # else
            #     dr = Δs[i+1][1]
            # end
            dl = i == 1 ? nothing : Δs[i-1][end]
            dr = i == n ? nothing : Δs[i+1][1]

            Δs[i] = partition(ds[i], _vs[i], dl, dr)
        end
    end
    @debug [length.(d) for d = Δs]

    r = map(zip(xs, Δs)) do (r, v)
        start = r[1][1]
        [start, (start + cumsum(reduce(vcat, v)))...]
    end
end
Base.getindex(::Nothing, i...) = nothing
Base.lastindex(::Nothing, i...) = nothing