using VideoIO
using ChainRulesCore

size2(a) = size(a, 1), size(a, 2)
size3(a) = size(a, 1), size(a, 2), size(a, 3)

function stuff(as)
    n = length(as)
    d = maximum(ndims.(as))
    L = l, w, h = maximum.([getindex.(size3.(as), i) for i = 1:3])
    i = j = k = round(n^(1 / d))
    if d == 2
        k = 1
    end
    if i * j * k < n
        i += 1
    end
    if i * j * k < n
        j += 1
    end
    if i * j * k < n
        k += 1
    end
    sz = (i, j, k)

    r = zeros(eltype(as[1]), sz .* (l, w, h))
    origins = zeros(Int, 3, n)
    for (li, (ix, a)) in enumerate(zip(CartesianIndices(sz), as))
        (i, j, k) = Tuple(ix)
        i = 1 + (i - 1) * l
        j = 1 + (j - 1) * w
        k = 1 + (k - 1) * h

        for (q, _l) = zip(ndims(a)+1:d, (l, w, h))
            # @show q, d
            a = stack(fill(a, L[q]))
        end
        ix = (i, j, k)
        r[range.(ix, ix + size3(a) - 1)...] = a
        origins[:, li] .= ix
    end
    r, origins
end
function unstuff(a, szs, origins)
    [
        begin
            m, n, o = [sz..., 1, 1]
            reshape(a[i:i+m-1, j:j+n-1, k:k+o-1], sz)
        end for ((i, j, k), sz) in zip(eachcol(origins), szs)
    ]
end

really(a::AbstractArray{<:Complex}) = vcat(real(a), imag(a))
really(a) = a

function _adjoint_reduce(f, ts, y, ulims)
    # encoder_options = (crf=23, preset="medium")
    encoder_options = (;)
    framerate = 24

    us, p, = y
    # _ass = []
    as = reduce(vcat, map(leaves.(us)) do u
        really.(u)
    end)

    A, origins = stuff(as)
    dummy = zeros(UInt8, size(A))
    depth = size(A, 3)
    sz0 = size(A)
    szs = size.(as)

    file = "video$(ts[1]).mp4"
    open_video_out(file, dummy[:, :, 1]; framerate, encoder_options) do writer
        for t = ts
            us, p, = y
            as = reduce(vcat, map(leaves.(us)) do u
                really.(u)
            end)
            # push!(_ass, reduce(vcat, leaves.(us)))

            A, = stuff(map(as, ulims) do a, (α, β)
                a = round.((a - α) / (β - α) * 255)
                a = max.(0, a)
                a = min.(255, a)
                convert.(UInt8, a)
            end)
            # @show size(A)
            for s = eachslice(A, dims=3)
                write(writer, s[:, :, 1])
            end
            y = f(y, t)
        end
    end
    return y, szs, origins, depth, file, sz0#, _ass
    #     open_video_out("video.mp4", stitch(stuff.(ass, szs)); framerate, encoder_options) do writer
end
adjoint_reduce(a...) = _adjoint_reduce(a...)[1]

function ChainRulesCore.rrule(config::RuleConfig{>:HasReverseMode}, ::typeof(adjoint_reduce), f, ts, y, ulims,)
    args = y[3:end]
    T = eltype(ts)
    us, p = y
    y, szs, origins, depth, file, sz0 = _adjoint_reduce(f, ts, y, ulims,)

    function adjoint_reduce_pullback(ūs)
        t̄s = zeros(T, length(ts))
        movie = VideoIO.openvideo(file, target_format=VideoIO.AV_PIX_FMT_GRAY8)
        virgin = true
        for (i, t) = zip(reverse(1:length(ts)), reverse(ts))
            vr = if i == 1
                VideoIO.seekstart(movie)
            else
                VideoIO.seek(movie, ((i - 1) * depth) / 24 - 0.001)
            end
            # @show VideoIO.counttotalframes(v)
            A = stack([reinterpret(UInt8, read(vr)) for (_,) = zip(1:depth,)])

            # global _szs, _file = szs, file
            as = map(unstuff(A, szs, origins), ulims) do a, (α, β)
                α + (β - α) / 255 * a |> T
            end
            # as = _ass[i]

            # T = typeof(as[1])
            # d = OrderedDict{Symbol,OrderedDict{Symbol,T}}()
            ks = keys(us[1])
            d = NamedTuple(Pair.(ks, as[1:length(ks)]))
            _y = ((d,), p, args...)

            j = 1 + length(ks)
            if (length(us)) == 2
                v = [[[[
                    begin
                        a = as[j]
                        j += 1
                        complex.(a[1:end÷2, :], a[end÷2+1:end, :])
                    end
                    for v = v] for v = v] for v = v] for v = us[2]]
                @assert j - 1 == sum(length.(leaves.(us)))
                _y = ((d, v), p, args...)
            end

            _, pb = rrule_via_ad(config, f, _y, t)
            # if virgin
            #     global _sfd = ȳ
            # end
            ȳ = pb(ūs)
            ūs, t̄ = ȳ[2:3]
            t̄s[i] = t̄
            # virgin && (global _sfd1 = y)
            # ȳ = (ȳ[1], ȳ[2], NoTangent(), NoTangent())
            virgin = false
            # ū, dp̄ = pb(ȳ)
            # p̄ += dp̄
        end
        # global _sfd2 = ȳ
        # println("p̄: ", p̄)
        # println("ū: ", ū)
        close(movie)
        return NoTangent(), NoTangent(), t̄s, ūs, NoTangent()
    end
    return y, adjoint_reduce_pullback
end