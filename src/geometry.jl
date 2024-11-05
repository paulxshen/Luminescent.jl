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
function tensorinv(a, lims, spacings)
    N = ndims(a)
    T = typeof(a)
    F = eltype(a)
    spacings = _downvec.(spacings, size(a))
    margin = [spacings[i][1] for i = 1:N]
    _a = pad(a, :replicate, margin)
    # A = cpu(A)

    v = [[
        begin
            li, ri = eachcol(lims[i])
            lj, rj = eachcol(lims[j])
            l = min.(li, lj) - 0.5
            _l = max.(li, lj) + 0.5
            Δ = _l - l
            start = round((l + 1) * margin) + 1
            global aa = [start, _l, l, spacings]
            downsample_by_range(_a[range.(start, start + sum.(spacings) + int((Δ - 1) * last.(spacings)) - 1)...], [[range(cum - space + 1, int(Δi .* space) + cum - space) for (cum, space) = zip(cumsum(spacing), spacing)] for (Δi, spacing) = zip(Δ, spacings)]) do a
                # n = ignore_derivatives() do
                n = if maximum(a) == minimum(a)
                    zeros(F, N)
                else
                    n = mean.([diff(a; dims) for dims in 1:N])
                    Z = norm(n)
                    if Z != 0
                        # if Z > 1e-2
                        n / Z
                    else
                        zeros(F, N)
                    end
                end
                # end
                P = n * n'
                P * mean(1 ./ a) + (I - P) / mean(a)
            end
        end for j = 1:i
    ] for i = 1:N]
    T.([
        begin
            if j > i
                j, i = i, j
            end
            getindex.(v[i][j], i, j)
        end for i = 1:N, j = 1:N
    ])
end