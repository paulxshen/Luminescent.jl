function partition(d, nres, v)
    n = ceil(Int, d * nres * sqrt(v))
    d / n, n
end
function transit(a, b)
    n = ceil(Int, log(1.5, b / a))
    r = (b / a)^(1 / n)
    a * r .^ (0:n)
end

function makemesh(mvs, bbox, nres)
    N = size(bbox, 1)
    rulers = [[(a, last(mvs)[2]), (b, nothing)] for (a, b) = eachrow(bbox)]
    for (m, v) = reverse(mvs)[2:end]
        vs = vertices(boundarybox(m))
        a = vs[1]
        b = vs[end]
        for (ruler, a, b) = zip(rulers, a, b)
            i = searchsortedfirst(ruler, a; by=first)
            j = searchsortedfirst(ruler, b; by=first)
            v1 = ruler[j-1][2]
            deleteat!(ruler, i:j-1)
            insert!(ruler, i, (a, v))
            insert!(ruler, j, (b, v1))
        end
    end
    deltas = map(rulers) do ruler
        map(enumerate(ruler[1:end-1])) do (i, (a, v))
            b = ruler[i+1][1]
            d = b - a

            v0 = v1 = d0 = d1 = nothing
            if i > 1
                d0 = a - ruler[i-1][1]
                v0 = ruler[i-1][2]
            end

            if i < length(ruler) - 1
                d1 = ruler[i+2][1] - b
                v1 = ruler[i+1][2]
            end

            l = i == 1 || v >= v0
            r = i == length(ruler) - 1 || v >= v1

            if l && r
                Δ, n = partition(d, nres, v)
                fill(Δ, n)
            elseif !l && !r
                Δ0, n0 = partition(d0, nres, v0)
                Δ1, n1 = partition(d1, nres, v1)
                L = transit(Δ0, Δ)
                R = transit(Δ1, Δ)
                s = d - (sum(L) + sum(R))
                if s > 0
                    Δs, ns = partition(s, nres, v)
                    vcat(L, fill(Δs, ns), reverse(R))
                else
                    Δ, n = partition(d, nres, max(v0, v1))
                    fill(Δ, n)
                end
            elseif l
                Δ1, n1 = partition(d1, nres, v1)
                R = transit(Δ1, Δ)
                s = d - sum(R)
                if s > 0
                    Δs, ns = partition(s, nres, v)
                    vcat(fill(Δs, ns), reverse(R))
                else
                    Δ, n = partition(d, nres, v1)
                    fill(Δ, n)
                end
            else
                Δ0, n0 = partition(d0, nres, v0)
                L = transit(Δ0, Δ)
                s = d - sum(L)
                if s > 0
                    Δs, ns = partition(s, nres, v)
                    vcat(L, fill(Δs, ns))
                else
                    Δ, n = partition(d, nres, v0)
                    fill(Δ, n)
                end
            end
        end
    end
    map(zip(rulers, deltas)) do (r, v)
        start = r[1]
        [start, (start + cumsum(vcat(reduce, v)))...]
    end
end

function samplemesh(meshvals, rulers)
    map(centroids(rulers)) do point
        for (i, (m, v)) = enumerate(meshvals)
            if isnothing(m) || sideof(point, m) == IN
                return v
            end
        end
    end
end