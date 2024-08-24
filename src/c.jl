using VideoIO
using ChainRulesCore
function getsize(n)
    i = j = round(sqrt(n))
    if i * j < n
        i += 1
    end
    if i * j < n
        j += 1
    end
    return (i, j)
end

size2(a) = size(a, 1), size(a, 2)
size3(a) = size(a, 1), size(a, 2), size(a, 3)

function stuff(as, sz=nothing, origins=nothing)
    n = length(as)
    l, w, h = maximum.([getindex.(size3.(as), i) for i = 1:3])
    i = j = k = round(n^(1 / 3))
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
        d, e, f = size3(a)
        r[i:i+d-1, j:j+e-1, k:k+f-1] = a
        origins[:, li] = [i, j, k]
    end
    r, origins
end
function unstuff(a, szs, origins)
    # m, n = size(a) ./ sz
    # map(CartesianIndices(sz)[1:k]) do (i, j)
    #     i = 1 + (i - 1) * m
    #     j = 1 + (j - 1) * n
    #     a[i:i+m-1, j:j+n-1]
    # end
    [
        begin
            m, n, o = [sz..., 1, 1]
            reshape(a[i:i+m-1, j:j+n-1, k:k+o-1], sz)
        end for ((i, j, k), sz) in zip(eachcol(origins), szs)
    ]
end
# function stitch(as)
#     r = zeros(eltype(as[1]), (sum(size.(as, 1)), maximum(size.(as, 2))))
#     i = 1
#     for a in as
#         m, n = size(a)
#         r[i:i+m-1, 1:n] = a
#         i += m
#     end
#     r
# end
# function unstitch(a, szs)
#     [a[i-m+1:i, 1:n] for (i, (m, n)) in zip(cumsum(size.(szs, 1)), szs)]
# end
complicateds = [false, true]
really(a::AbstractArray{<:Complex}) = vcat(real(a), imag(a))
really(a) = a
function _adjoint_reduce(f, ts, y, ls)
    encoder_options = (crf=23, preset="medium")
    framerate = 24

    us, p, = y
    as = reduce(vcat, map(leaves.(us)) do u
        really.(u)
    end)

    szs = size.(as)
    A, origins = stuff(as)
    dummy = zeros(UInt8, size(A))
    depth = size(A, 3)
    sz = size(A)

    file = "video$(round(length(ts))).mp4"
    open_video_out(file, dummy[:, :, 1]; framerate, encoder_options) do writer
        for t = ts
            us, p, = y
            as = reduce(vcat, map(leaves.(us)) do u
                really.(u)
            end)

            A, = stuff(map(as, ls) do a, (α, β)
                    # @show a, α, β
                    if α == β
                        α -= 1
                        β += 1
                    end
                    if eltype(a) <: Complex
                        a = vcat(reim(a)...)
                    end
                    a = round.((a - α) / (β - α) * 255)
                    a = max.(0, a)
                    a = min.(255, a)
                    convert.(UInt8, a)
                end, sz, origins)
            for s = eachslice(A, dims=3)
                write(writer, s)
            end
            y = f(y, t)
        end
    end
    return y, szs, origins, depth, file
    #     open_video_out("video.mp4", stitch(stuff.(ass, szs)); framerate, encoder_options) do writer
end
adjoint_reduce(a...) = _adjoint_reduce(a...)[1]

function ChainRulesCore.rrule(config::RuleConfig{>:HasReverseMode}, ::typeof(adjoint_reduce), f, ts, y, ls,)
    args = y[3:end]
    us, p = y
    y, szs, origins, depth, file = _adjoint_reduce(f, ts, y, ls,)
    # ns = length.(leaves.(u))
    # grid_sizes = getsize.(ns)
    # block_sizes = [gsz .*
    #                begin
    #     v = size2.(as)
    #     ms = getindex.(v, 1)
    #     ns = getindex.(v, 2)
    #     maximum(ms), maximum(ns)
    # end
    #                for (gsz, as) = zip(grid_sizes, leaves.(u))]
    # nds = [ndims(leaves(u)[1]) for u = u]

    function adjoint_reduce_pullback(ȳ)
        # ȳ = ȳ[2]
        movie = VideoIO.openvideo(file)
        # for (img, t) in zip(movie, reverse(ts))
        # ūs = ȳ[1:end-1]
        # p̄ = ȳ[end]
        # ū = ȳ
        # y = collect(ȳ)
        # if y[1] != ZeroTangent()
        #     y[1] = NamedTuple(pairs(y[1]))
        # end
        # ū = Tuple(y)
        virgin = true
        for (i, t) = zip(reverse(1:length(ts)), reverse(ts))
            v = if i == 1
                # error()
                VideoIO.seekstart(movie)
            else
                VideoIO.seek(movie, ((i - 1) * depth) / 24 - 0.001)
            end
            # @show i, t #length(v)
            A = stack([reinterpret(UInt8, a) for (_, a) = zip(1:depth, v)])

            as = map(unstuff(A, szs, origins), ls) do a, (α, β)
                α + (β - α) / 255 * a
            end

            # T = typeof(as[1])
            # d = OrderedDict{Symbol,OrderedDict{Symbol,T}}()
            ks = keys(us[1])
            d = NamedTuple(Pair.(ks, as[1:length(ks)]))
            _y = ((d,), p, args...)

            i = 1 + length(ks)
            if (length(us)) == 2
                v = [[[[
                    begin
                        a = as[i]
                        i += 1
                        complex.(a[1:end÷2, :], a[end÷2+1:end, :])
                    end
                    for v = v] for v = v] for v = v] for v = us[2]]
                _y = ((d, v), p, args...)
            end

            _, pb = rrule_via_ad(config, f, _y, t)
            # if virgin
            #     global _sfd = ȳ
            # end
            ȳ = pb(ȳ)
            ȳ = ȳ[2]
            # virgin && (global _sfd1 = y)
            ȳ = (ȳ[1], ȳ[2], NoTangent(), NoTangent())
            virgin = false
            # ū, dp̄ = pb(ȳ)
            # p̄ += dp̄
        end
        # global _sfd2 = ȳ
        # println("p̄: ", p̄)
        # println("ū: ", ū)
        close(movie)
        return NoTangent(), NoTangent(), NoTangent(), ȳ, NoTangent()
    end
    return y, adjoint_reduce_pullback
end