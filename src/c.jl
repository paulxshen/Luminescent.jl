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
function adjoint_reduce(f, ts, u, p, ls, withmovie=false)
    encoder_options = (crf=23, preset="medium")
    framerate = 24

    as = reduce(vcat, [[
        begin
            if eltype(a) <: Complex
                a = vcat(real(a), imag(a))
            end
            a
        end for a = v
    ] for v = leaves.(u)])
    szs = size.(as)
    A, origins = stuff(as)
    dummy = zeros(UInt8, size(A))
    depth = size(A, 3)
    sz = size(A)

    file = "video$(round(length(ts))).mp4"
    open_video_out(file, dummy[:, :, 1]; framerate, encoder_options) do writer
        for t = ts
            if withmovie
                A, = stuff(map(reduce(vcat, leaves.(u)), ls) do a, (α, β)
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
            end
            u = f(u, p, t)
        end
    end
    if withmovie
        return u, szs, origins, depth, file
    end
    #     open_video_out("video.mp4", stitch(stuff.(ass, szs)); framerate, encoder_options) do writer
    u
end

function ChainRulesCore.rrule(config::RuleConfig{>:HasReverseMode}, ::typeof(adjoint_reduce), f, ts, u, p, ls, with)
    u0 = deepcopy(u)
    y, szs, origins, depth, file = adjoint_reduce(f, ts, u, p, ls, true)
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
        movie = VideoIO.openvideo(file)
        # for (img, t) in zip(movie, reverse(ts))
        p̄ = 0
        # ū = ȳ
        y = collect(ȳ)
        if y[1] != ZeroTangent()
            y[1] = NamedTuple(pairs(y[1]))
        end
        ū = Tuple(y)

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

            i = 1
            T = typeof(as[1])
            d = OrderedDict{Symbol,OrderedDict{Symbol,T}}()
            u1 = 0
            for k1 = keys(u0[1])
                d[k1] = OrderedDict{Symbol,T}()
                for k2 = keys(u0[1][k1])
                    d[k1][k2] = as[i]
                    i += 1
                end
                u1 = (d,)
            end

            if (length(u0)) == 2
                v = [[[[
                    begin
                        a = as[i]
                        i += 1
                        complex.(a[1:end÷2, :], a[end÷2+1:end, :])
                    end
                    for v = v] for v = v] for v = v] for v = u0[2]]
                u1 = (d, v,)
            end

            _, pb = rrule_via_ad(config, f, u1, p, t)
            ū, dp̄ = pb(ū)[2:3]
            p̄ += dp̄
        end
        println("p̄: ", p̄)
        println("ū: ", ū)
        close(movie)
        return NoTangent(), NoTangent(), NoTangent(), ū, p̄, NoTangent(), NoTangent()
    end
    return y, adjoint_reduce_pullback
end