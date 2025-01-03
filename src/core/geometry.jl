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
function tensorinv(a::T, lims, spacings) where {T}
    N = ndims(a)
    F = eltype(a)


    spacings = _downvec.(spacings, size(a))
    margin = [spacings[i][1] for i = 1:N]
    _a = pad(a, :replicate, margin)
    v = [[
        begin
            li, ri = eachcol(lims[i])
            lj, rj = eachcol(lims[j])
            l = min.(li, lj) - 0.5
            _l = max.(li, lj) + 0.5
            Δ = _l - l

            start = round((l + 1) * margin) + 1
            stop = start + sum.(spacings) + int((Δ - 1) * last.(spacings)) - 1
            __a = _a[range.(start, stop)...]
            Is = [[range(cum - space + 1, int(Δi .* space) + cum - space) for (cum, space) = zip(cumsum(spacing), spacing)] for (Δi, spacing) = zip(Δ, spacings)]


            As = downsample_by_range(identity, __a, Is)
            M = mean.(As)
            invM = map(As) do a
                mean(1 ./ a)
            end
            V = [sum.(diff.(As; dims)) for dims in 1:N]
            Z = sqrt.(sum([a .^ 2 for a in V])) + F(0.001)
            V = V ⊘ [Z]

            Pij = V[i] .* V[j]
            Pij .* invM + ((i == j) - Pij) ./ M
        end for j = 1:i
    ] for i = 1:N]

    T.([j <= i ? v[i][j] : v[j][i] for i = 1:N, j = 1:N])
end