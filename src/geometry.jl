function apply_subpixel_averaging(geometry, fieldlims)
    namedtuple([k => begin
        lrs = ignore_derivatives() do
            E = r"E.*"
            H = r"H.*"
            f = (; ϵ=E, σ=E, μ=H, m=H, invϵ=E)[k]
            fieldlims(f)
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
                pad(a, geometry_padvals[k][1], eachcol(geometry_padamts[k] * ratio)...)
            else
                [a]
            end
        end for k = keys(geometry)])
end

function tensorinv(a, fieldlims, ratio,)
    @assert ratio > 1
    N = ndims(a)
    T = typeof(a)
    F = eltype(a)
    sz = size(a) .÷ ratio
    margin = ratio ÷ 2
    A = pad(a, :replicate, margin)
    A = cpu(A)

    v = [
        begin
            lr = cpu(lr)
            l, r = eachcol(lr)
            start = Int.((l - 0.5) * ratio + margin) + 1
            finish = size(A) + Int.((r - sz + 0.5) * ratio - margin)
            # @show start, finish, size(A), margin, l, r, sz
            _a = A[range.(start, finish)...]
            # _a = crop(A, (l - 0.5) * ratio + margin, margin - (r - sz + 0.5) * ratio)
            downsample(_a, ratio) do a
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
        end for lr = fieldlims(r"E.*")
    ]
    T.([getindex.(v[j], i, j) for i = 1:N, j = 1:N])
    # @assert !any(iszero, a)
    # @assert !any(isnan, a)
    # # return [i == j ? 1 ./ a : zeros(F, size(a)) for i = 1:N, j = 1:N]
    # [i == j ? 1 ./ downsample(a, ratio) : zeros(F, size(a) .÷ ratio) for i = 1:N, j = 1:N]
end