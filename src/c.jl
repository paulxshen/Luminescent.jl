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

function stuff(as, sz)
    m, n = size(as[1])
    r = zeros(eltype(as[1]), sz)
    for ((i, j), a) in zip(CartesianIndices(sz), as)
        m, n = size(as[1])
        i = 1 + (i - 1) * m
        j = 1 + (j - 1) * n
        r[i:i+m-1, j:j+n-1] = a
    end
    r
end
function unstuff(a, sz, k)
    m, n = size(a)
    map(CartesianIndices(sz)[1:k]) do (i, j)
        i = 1 + (i - 1) * m
        j = 1 + (j - 1) * n
        a[i:i+m-1, j:j+n-1]
    end
end

complicateds = [false, true]
function adjoint_reduce(f, ts, u, lss; withmovie=false)
    encoder_options = (crf=23, preset="medium")
    framerate = 24
    ass = [[
        begin
            if eltype(a) <: Complex
                a = vcat(real(a), imag(a))
            end
            a
        end for a = v
    ] for v = leaves.(u)]
    szs = getsize.(length.(ass))
    open_video_out("video.mp4", stitch(stuff.(ass, szs)); framerate, encoder_options) do writer
        for t = ts
            u = f(u, t)
            withmovie && write(writer, stitch([stuff(map(leaves(u), ls) do a, (p, q)
                    a = (a - p) / (q - p) * 256
                    a = max.(0, a)
                    a = min.(255, a)
                    round.(UInt8, a)
                end, sz) for (u, ls, sz) in zip(u, lss, szs)]))
        end
    end

    if withmovie
        return u, VideoIO.openvideo("video.mp4")
    end
    u
end

function ChainRulesCore.rrule(config::RuleConfig{>:HasReverseMode}, ::typeof(adjoint_reduce), f, ts, u, lss)
    y, m = adjoint_reduce(f, t, u, lss; withmovie=true)
    ns = length.(leaves.(u))
    szs = getsize.(ns)
    function adjoint_reduce_pullback(ȳ)
        f̄ = NoTangent()
        t̄s = NoTangent()

        for (img, t) in zip(reverse(m), reverse(ts))
            _u = [
                begin
                    as = map(unstuff(a, sz, n), ls) do a, (p, q)
                        r = p + (q - p) / 256 * a
                        if complicated
                            r = r[1:end÷2, :] + r[end÷2+1:end, :] * im
                        end
                        r
                    end
                    i = 1
                    d = OrderedDict()
                    for k1 = keys(u)
                        d[k1] = OrderedDict()
                        for k2 = keys(u[k1])
                            d[k1][k2] = as[i]
                            i += 1
                        end
                    end
                    d
                end for (ls, n, sz, a, complicated) = zip(lss, ns, szs, unstitch(img, getindex.(szs, 1)), complicateds)
            ]
            _, pb = rrule_via_ad(config, f, _u, t)
            ȳ, _ = pb(ȳ)
        end
        return f̄, t̄s, ȳ, NoTangent()
    end
    return y, adjoint_reduce_pullback
end