
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
function tensorinv(a, pec, field_lims, spacings, prop=:ϵ)
    if !any(pec) && isa(a, Number)
        return 1 / a
    end

    T = eltype(a)
    N = ndims(pec)
    sz = size(pec)
    if prop == :ϵ
        lims = values(field_lims(r"E.*"))
    elseif prop == :μ
        lims = values(field_lims(r"H.*"))
    end

    spacings = _downvec.(spacings, sz)
    # global _t = spacings, lims

    ratio = minimum(minimum.(spacings))
    if ratio >= 4
        border = 0
    else
        border = 1
    end
    @show border

    margin = map(spacings) do s
        (1 + border) * max(first(s), last(s))
    end
    ap = isa(a, AbstractArray) ? pad(a, :replicate, margin) : a
    pecp = isa(pec, AbstractArray) ? pad(pec, :replicate, margin) : pec

    v = map(1:N) do i
        map(1:i) do j
            li, ri = eachcol(lims[i])
            lj, rj = eachcol(lims[j])
            l = min.(li, lj) - (0.5 + border)
            _l = max.(li, lj) + 0.5 + border
            Δ = _l - l
            start = round((l + 1 + border) * margin) + 1

            # global _d = Δ, spacings, lims
            I = range.(start, start + sum.(spacings) + int((Δ - 1) * last.(spacings)) - 1)
            api = ap(I...)
            pecpi = pecp(I...)


            ranges = [[
                range(cum - space + 1, int(Δi .* space) + cum - space)
                for (cum, space) = zip(cumsum(spacing), spacing)
            ] for (Δi, spacing) = zip(Δ, spacings)]

            downsample_by_range((api, pecpi), ranges) do a, pec
                # amin, amax = extrema(a)
                # if all(pec)
                #     return 0 |> T
                # elseif amin == amax && !any(pec)
                #     return (i == j) / amax |> T
                # else
                begin
                    v = a
                    if prop == :μ
                        if any(pec)
                            v = pec
                        end
                    end

                    if isa(v, Number)
                        Pij = 0
                    else
                        n = imnormal(v)
                        Pij = n[i] * n[j]
                    end

                    # if any(pec)
                    #     f = sum(pec) / prod(size(pec))
                    #     if prop == :ϵ
                    #         return Pij / amin * f |> T
                    #         # return Pij / amin |> T
                    #     elseif prop == :μ
                    #         return ((i == j) - Pij) / amin * f |> T
                    #         # return ((i == j) - Pij) / amin |> T
                    #     else
                    #         error("prop must be either ϵ or μ")
                    #     end
                    # else
                    if prop == :ϵ
                        Pij * mean(1 ./ a) + ((i == j) - Pij) / mean(a)
                    elseif prop == :μ
                        if any(pec)
                            ((i == j) - Pij)
                        else
                            (i == j)
                        end
                    else
                        error("prop must be either ϵ or μ")
                    end
                end |> T
            end
            # end
        end
    end

    [j <= i ? v[i][j] : v[j][i] for i = 1:N, j = 1:N]
end



function pecfy(a::AbstractArray{T}, field_lims, spacings) where {T}
    N = ndims(a)
    lims = values(field_lims(r"E.*"))
    spacings = _downvec.(spacings, size(a))
    rulers = cumsum.(spacings)

    pec = a .>= (PECVAL - TOL)
    d = DefaultDict{T,BitArray{N}}(passkey=true) do k
        @show k
        zeros(Bool, size(pec))
    end
    # d = Dict()
    # for _ = 1:1
    for (i, lims) = enumerate(lims)
        I = map(eachrow(lims), rulers) do (l, r), ruler
            map(range(l, r, length=int(r - l + 1))) do i
                start = 1 + getindexs(ruler, i) |> int
                stop = getindexs(ruler, i + 1) |> int

                start:stop
            end
        end

        sel = (i .== 1:N)
        ignore_derivatives() do
            [
                begin
                    _pec = pec[I...]
                    p, q = extrema(_pec)
                    if p != q
                        m = mean(_pec)
                        kmin, kmax = extrema(a[I...])
                        kmode = m > 0.5 ? kmax : kmin
                        k = nothing
                        if m < 0.2
                            k = kmin
                        elseif m > 1 - 0.2
                            k = kmax
                            # !haskey(d, k) && d[k] = zeros(Bool, size(pec))
                        else
                            n = imnormal(a)
                            if all(iszero, n) || abs(n ⋅ sel) > cosd(45)
                                k = kmode
                            end
                        end
                        if !isnothing(k)
                            # d[k][I...] .= true
                            pec[I...] .= true
                        end
                    end
                end for I = Base.product(I...)
            ]
        end
    end
    # end

    a = a .* ((d |> values |> sum) .== 0) + sum(pairs(d)) do (k, v)
        k * v
    end
    a, pec
end