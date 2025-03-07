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
function pad_geometry(geometry, geometry_padvals, geometry_padamts, ratio=1)
    namedtuple([
        k => begin
            a = geometry[k]
            if isa(a, AbstractArray) && k in keys(geometry_padvals) && k in keys(geometry_padamts)
                pad(a, geometry_padvals[k][1], eachcol(geometry_padamts[k] .* ratio)...)
            else
                a
            end
        end for k = keys(geometry)])
end

_size(s::Real, n) = int(s * n)
_size(s, _) = int(sum(s))

function tensorinv(a::AbstractArray{T}, field_lims, spacings) where {T}
    N = ndims(a)
    # lims = values(field_lims(r"E.*"))
    lims = @ignore_derivatives [field_lims("E$v") for v = collect("xyz")[1:N]]
    spacings = _downvec.(spacings, size(a))

    border = 0
    margin = map(spacings) do s
        (1 + border) * max(first(s), last(s))
    end
    ap = pad(a, :replicate, margin)

    v = map(1:N) do i
        map(1:i) do j
            li, ri = eachcol(lims[i])
            lj, rj = eachcol(lims[j])
            l = min.(li, lj) - 0.5
            _l = max.(li, lj) + 0.5
            Δ = _l - l
            start = round((l + 1) * margin) + 1

            # global _d = Δ, spacings, lims

            api = ap[range.(start, start + sum.(spacings) + int((Δ - 1) * last.(spacings)) - 1)...]
            ranges = [[
                range(cum - space + 1, int(Δi .* space) + cum - space)
                for (cum, space) = zip(cumsum(spacing), spacing)
            ] for (Δi, spacing) = zip(Δ, spacings)]

            downsample_by_range(api, ranges) do a
                p, q = extrema(a)
                p == q & return (i == j) / p

                n = imnormal(v)
                Pij = n[i] * n[j]
                Pij * mean(1 ./ a) + ((i == j) - Pij) / mean(a)
            end |> T
        end
    end

    [j <= i ? v[i][j] : v[j][i] for i = 1:N, j = 1:N]
end

function tensorinv(meshvals::AbstractVector{Tuple}, rulers)
    N = length(rulers)
    ns = length.(rulers)
    a = map(Base.product(range.(OneTo.(ns - 1))...)) do I
        start = getindex.(rulers, I)
        stop = getindex.(rulers, I .+ 1)
        Δ = stop - start

        start -= Δ / 2
        stop += Δ / 2
        start .= max.(start, first.(rulers))
        stop .= min.(stop, last.(rulers))

        box = Box(start, stop)
        point = centroid(box)
        hits = fill(false, length(meshvals))
        for (i, (m, v)) = enumerate(meshvals)
            if !isnothing(m) && intersects(box, m)
                hits[i] = true
            elseif isnothing(m) || sideof(point, m) == IN
                return v
            end
        end

        n = 4
        δ = Δ / n
        a = map(Base.product(range.(start + δ / 2, stop - δ / 2, n)...)) do p
            for (m, v) = meshvals[hits]
                (isnothing(m) || sideof(p, m) != OUT) && return v
            end
            default
        end

        n = imnormal(a)
        P = n * n'
        P * mean(1 ./ a) + (LinearAlgebra.I - P) / mean(a)
    end

    v = map(1:N) do i
        map(1:i) do j
            map(a) do v
                v(i, j)
            end
        end
    end
    [j <= i ? v[i][j] : v[j][i] for i = 1:N, j = 1:N]
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