function apply_subpixel_averaging(geometry, field_lims)
    namedtuple([k => begin
        lrs = ignore_derivatives() do
            E = r"E.*"
            H = r"H.*"
            f = (; ϵ=E, σ=E, μ=H, m=H, invϵ=E)[k]
            field_lims(f)
        end
        @assert length(lrs) > 0
        a = geometry[k]
        map(lrs) do lr
            crop(a, lr[:, 1] - 0.5, size(a) - 0.5 - lr[:, 2])
        end
        # end for (k, a) = pairs(geometry)])
    end for k = keys(geometry)])
end

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

function tensorinv(a, field_lims, field_spacings)
    N = ndims(a)
    T = typeof(a)
    F = eltype(a)
    margin = round([field_spacings[i](1)[1] for i = 1:N] / 2)
    _a = pad(a, :replicate, margin, 0)
    # A = cpu(A)

    v = [
        begin
            fs = round.(vec.(getindex.(field_spacings, k)))
            l, r = eachcol(field_lims[k])
            start = round(l .* margin) + 1
            sz = sum.(fs)
            global aa = [start, sz, fs, _a]
            downsample(_a[range.(start, start + sz - 1)...], fs) do a
                n = ignore_derivatives() do
                    if maximum(a) == minimum(a)
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
                end
                P = n * n'
                P * mean(1 ./ a) + (I - P) / mean(a)
            end
        end for k = keys(field_lims(r"E.*"))]
    T.([getindex.(v[j], i, j) for i = 1:N, j = 1:N])
end