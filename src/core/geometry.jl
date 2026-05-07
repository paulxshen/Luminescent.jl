function apply_subpixel_averaging(geometry::Number, args...)
    geometry
end
function apply_subpixel_averaging(geometry, field_lims)
    namedtuple([k => begin
        a = geometry[k]
        if isa(a, Number)
            a
        else
            E = r"E.*"
            H = r"H.*"
            f = (; ϵ=E, σ=E, μ=H, m=H, invϵ=E)[k]
            lrs = field_lims(f)
            @assert length(lrs) > 0
            map(lrs) do lr
                crop(a, lr[:, 1] - 0.5, size(a) - 0.5 - lr[:, 2])
            end
        end
        # end for (k, a) = pairs(geometry)])
    end for k = keys(geometry)])
end

pad_geometry(geometry::Number, args...) = geometry
function pad_geometry(geometry, geometry_padvals, padamts)
    namedtuple([
        k => begin
            v = geometry[k]
            if v isa Number || isnothing(v) || k ∉ keys(geometry_padvals) #|| k ∉ keys(geometry_padamts)
                v
            else
                padval = geometry_padvals[k]
                # padamt = geometry_padamts[k]
                if k == :Cn
                    map(v) do v
                        map(v) do v
                            map(v) do a
                                pad(a, eachcol(padval)..., eachcol(padamts)...)
                            end
                        end
                    end
                else
                    map(v) do a
                        pad(a, eachcol(padval)..., eachcol(padamts)...)
                    end
                end
            end
        end for k = keys(geometry)])
end

inveps(x) = isPEC(x) ? 0 : 1 / x

function supersamplegrid(a::AbstractArray{S,N}, vals, deltas, ratio, offsets, invtensor) where {S,N}
    isovals = mean.(diag.(vals))
    invvals = inveps.(vals)
    @debug vals, isovals, invvals
    if invtensor
        sz = int.(size(a) / ratio) - 2
        A = map(1:N) do i
            map(1:i) do j
                # offset = int.(ratio * (offsets(i) + offsets(j)) / 2)
                offset = @ignore_derivatives (offsets(i) + offsets(j)) / 2
                tmap(CartesianIndices(sz)) do I
                    I = Tuple(I)
                    start = int.((I + offset) * ratio) + 1
                    stop = start + ratio - 1
                    d = call.(deltas, I + 1)

                    E = a[(:).(start, stop)...]
                    # T = eltype(isovals)
                    # r = getindex.((isovals,), E)

                    # r = a[(:).(start, stop)...]
                    # p, q = @ignore_derivatives extrema(r)
                    # if isPEC(q) || isPEC(p)
                    #     (i != j) && return 0 |> T
                    #     v = r((size(r) / 2 + 0.5)...)
                    #     isPEC(v) && return 0 |> T
                    #     return 1 / v
                    # end
                    # nc

                    # return (i == j) / mean(r)
                    p, q = extrema(E)
                    p == q && return (i == j) * invvals[int(E[1] + 1)](i, j)

                    v0 = isovals[1](i, j)
                    v1 = isovals[2](i, j)
                    r = v0 + (v1 - v0) * E
                    n = imnormal(r, d)
                    Pij = n[i] * n[j]
                    Pij * mean(1 ./ r) + ((i == j) - Pij) / mean(r)
                end
            end
        end
        [j <= i ? A[i][j] : A[j][i] for i = 1:N, j = 1:N]
    else
        # return map(offsets) do offset
        sz = int.(size(a) / ratio)
        offsets = @ignore_derivatives collect(_values(offsets))
        return map(offsets) do offset
            tmap(CartesianIndices(sz)) do I
                I = Tuple(I)
                # I = @ignore_derivatives (I + offset) * ratio + (ratio + 1) / 2
                # a(I...)
                start = int.((I + offset - 1) * ratio) + 1
                stop = start + ratio - 1
                start = max.(start, 1)
                stop = min.(stop, size(a))
                E = a[(:).(start, stop)...]
                # r = getindex.((isovals,), E)
                v0 = isovals[1]
                v1 = isovals[2]
                r = v0 + (v1 - v0) * E
                mean(r)
            end
        end
    end
end

function vec3(p::Meshes.Point)
    p = coords(p)
    # paramdim(p) == 2 && return [p.x, p.y]
    # paramdim(p) == 3 && 
    return SVector{3}(ustrip([p.x |> FF, p.y |> FF, p.z |> FF]))
end
vec3(v) = SVector{3}(v)

struct Searcher
    tree
    ps
    ns
    bbox
    isbox
end
function Searcher(m, d, tol)
    x = boundingbox(m)
    a = ustrip.(getfield(coords(x.min), :coords)) |> vec3
    b = ustrip.(getfield(coords(x.max), :coords)) |> vec3
    bbox = [a b]

    cs = vec3.(centroid.(m))
    isbox = all(cs) do p
        any(zip(a, b, p)) do (a, b, p)
            isapprox(p, a; atol=tol) || isapprox(p, b; atol=tol)
        end
    end

    if isbox
        v = map(enumerate(eachrow(bbox))) do (i, (zmin, zmax))
            I = collect(1:3)
            deleteat!(I, i)
            xmin, xmax = bbox[I[1], :]
            ymin, ymax = bbox[I[2], :]
            nx = ceil(Int, (xmax - xmin) / d)
            ny = ceil(Int, (ymax - ymin) / d)
            dx = (xmax - xmin) / nx
            dy = (ymax - ymin) / ny
            v = reduce(hcat, collect.(Base.product(
                range(xmin + dx / 2, xmax - dx / 2, nx),
                range(ymin + dy / 2, ymax - dy / 2, ny))))
            push!(I, i)
            n = size(v, 2)

            ps = reduce(hcat, map((zmin, zmax)) do z
                w = vcat(v, fill(z, n)')
                stack(eachrow(w)[invperm(I)])'
            end)
            ns = reduce(hcat, map((-1, 1)) do v
                a = zeros(FF, 3, n)
                a[i, :] .= v
                a
            end)

            ps, ns
        end
        ps = reduce(hcat, first.(v))
        ns = reduce(hcat, last.(v))
    else
        sampler = MinDistanceSampling(d)
        tris = tmap(collect(m)) do f
            c = vec3(centroid(f))
            try

                ps1 = sample(f, sampler) .|> vec3
            catch
                ps1 = [c]
            end
            if isempty(ps1)
                ps1 = [c]
            end
            # v = reduce(hcat, vec3.(vertices(f)))
            # ps1 = [c, ps1...]
            n = vec3(ustrip(normal(f)))
            # p = perimeter(f) |> ustrip
            (; c, n, ps=ps1)
        end
        ps = getindex.(tris, :ps)
        ns = reduce(hcat, reduce(vcat, fill.(getindex.(tris, :n), length.(ps))))
        ps = reduce(hcat, reduce(vcat, ps))
    end

    tree = BallTree(ps)
    @debug "Searcher tol=$tol isbox=$isbox"
    Searcher(tree, ps, ns, bbox, isbox)
end


function Base.in(p, s::Searcher, tol)
    @unpack tree, ps, ns, bbox, isbox = s
    for (v, (a, b)) = zip(p, eachrow(bbox))
        (v < a - tol || v > b + tol) && return false
    end
    isbox && return true

    i, = nn(tree, p)
    P = ps[:, i]
    N = ns[:, i]
    (p - P) ⋅ N ≤ tol
end

function samplemesh(point, funnel, searchers, vals, tol)
    @assert length(searchers) + 1 == length(vals)
    bg = vals[end]
    I = findall(point, funnel)
    for (searcher, val) in zip(searchers[I], vals[I])
        # point ∈ searcher && return val
        Base.in(point, searcher, tol) && return val
    end
    bg
end

function proj(point, delta, funnel, searchers, bbox)
    vals = 1:length(searchers)+1
    I1 = findall(point - delta / 2, point + delta / 2, funnel)
    bg = vals[end]
    isempty(I1) && return nothing, nothing, nothing, bg, bg

    tol = mean(delta) / 100
    ds = ns = ls = []
    vin = nothing
    for (val, searcher) in zip(vals[I1], searchers[I1],)
        I = inrange(searcher.tree, point, maximum(delta / 2))
        ns = searcher.ns[:, I]
        ds = eachcol(point .- searcher.ps[:, I]) .⋅ eachcol(ns)
        ls = vec(abs.(delta' * ns))
        # ls = map(eachcol(ns)) do n
        #     norm(delta .* n)
        # end

        I = findall(abs.(ds) .<= ls / 2)
        ds = ds[I]
        ls = ls[I]
        ns = ns[:, I]

        for (i, d) = enumerate(ds)
            if d <= 0
                vin = val
                ds[i] *= -1
            else
                ns[:, i] *= -1
            end
        end
        if isnothing(vin)
            vin = samplemesh(point, funnel, searchers, vals, tol)
        end
        !isempty(ds) && break
    end

    if isempty(ds)
        v = samplemesh(point, funnel, searchers, vals, tol)
        return nothing, nothing, nothing, v, v
    end

    # if length(ds) > 1
    #     I = []
    #     for i = 1:length(ds)
    #         n = ns[:, i]
    #         new = true
    #         for _i = I
    #             _n = ns[:, _i]
    #             if all(<(0.01), abs.(n - _n))
    #                 new = false
    #                 break
    #             end
    #         end
    #         if new
    #             push!(I, i)
    #         end
    #     end
    #     ns = ns[:, I]
    #     ds = ds[I]
    #     ls = ls[I]
    # end

    # single = true
    # if length(ds) > 1
    #     ws = 1 - 2ds ./ ls
    #     @assert all(>(-0.01), ws)
    #     ws = max.(ws, 0)
    #     s = sum(ws)
    #     if s > 0
    #         single = false
    #         ws ./= s

    #         d = ds ⋅ ws
    #         l = ls ⋅ ws
    #         n = ns * ws
    #     end
    # end

    # if single
    #     d = ds[1]
    #     l = ls[1]
    #     n = ns[:, 1]
    # end

    i = 1
    if length(ds) > 1
        i = argmin(ds ./ ls)
    end
    d = ds[i]
    l = ls[i]
    n = ns[:, i]

    p = point + n * l / 2
    p += 2(max.(0, bbox[:, 1] - p) + min.(0, bbox[:, 2] - p))
    vout = samplemesh(p, funnel, searchers, vals, tol)

    if vin == vout
        d = n = l = nothing
    end
    d, n, l, vin, vout
end

function supersamplemesh(centers, deltas, funnel, searchers, bbox, vals, offsets, invtensor=false)
    N = length(centers)
    sz = Tuple(length.(centers))
    isovals = mean.(diag.(vals))
    invvals = inveps.(vals)
    if invtensor
        A = map(1:N) do i
            map(1:i) do j
                offset = @ignore_derivatives (offsets(i) + offsets(j)) / 2
                tmap(CartesianIndices(sz)) do I
                    I = Tuple(I)
                    origin = getindex.(centers, I)
                    delta = getindex.(deltas, I)
                    origin += offset .* delta
                    d, n, l, ein, eout = proj(origin, delta, funnel, searchers, bbox)
                    isnothing(d) && return (i == j) * invvals[ein](i, j)

                    vin = isovals[ein]
                    vout = isovals[eout]
                    isPEC(vin) && return 0
                    isPEC(vout) && return (i == j) * invvals[ein](i, j)
                    Pij = n[i] * n[j]
                    win = 0.5 + d / l |> FF
                    @assert 0.49 <= win <= 1.01 win
                    wout = 1 - win
                    Pij * (win / vin + wout / vout) + ((i == j) - Pij) / (win * vin + wout * vout)
                end
            end
        end

        return [j <= i ? A[i][j] : A[j][i] for i = 1:N, j = 1:N]
    else
        map(offsets) do offset
            tmap(CartesianIndices(sz)) do I
                I = Tuple(I)
                origin = getindex.(centers, I)
                delta = getindex.(deltas, I)
                origin += offset .* delta
                tol = mean(delta) / 100
                samplemesh(origin, funnel, searchers, vals, tol)
            end
        end
    end
end



function downsamplefield(a::AbstractArray{T}, field_lims, spacings) where {T}
    N = ndims(a)
    lims = values(field_lims(r"E.*"))
    spacings = _downvec.(spacings, size(a))
    rulers = cumsum.(spacings)

    map(lims) do lims
        I = map(eachrow(lims), rulers) do (l, r), ruler
            map(range(l, r, length=int(r - l + 1))) do i
                start = 1 + getindexs(ruler, i) |> int
                stop = getindexs(ruler, i + 1) |> int

                start:stop
            end
        end
        downsample_by_range(mean, a, I)
    end
end