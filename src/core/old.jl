function _supersamplemesh(T, N, searchers, vals, bg, offset, centers, Δs, _inbound, tensor, z, i=0, j=0)
    @unpack tree, tris, tri2body = searchers
    tmap(collect(zip(Base.product(centers...), Base.product(Δs...)))) do (origin, Δ)
        origin += offset .* Δ
        if N == 2
            @assert !isnothing(z)
            origin = [origin[1], origin[2], z]
            Δ = [Δ[1], Δ[2], (Δ[1] + Δ[2]) / 2]
        end

        vedge = edge_wt = center_normal = nothing
        R = (0.65)
        R *= Δ

        tol = minimum(Δ) / 100
        is = inrange(tree, origin, maximum(R))
        # v = map(unique(getindex.((p2t,), is))) do i
        v = map(is) do i
            c, n, v = tris[i]
            v = v .- origin
            !raytri(n, v) && return nothing
            d = -n ⋅ (c - origin)
            r = norm(d * n ./ R)
            r > 1 && return nothing
            if d ≈ 0
                # d = sign(findfirst(x -> !(x ≈ 0), n)) * tol
                d = -tol
            end
            d, r, n
        end |> filter(!isnothing)

        if !isempty(v)
            v = map(v) do (d, r, n)
                Z = norm(n[1:N])
                Z == 0 && return nothing

                h = 1 - r
                if N == 3
                    w = h^2 * (3 - h) / 4
                else
                    w = (acos(1 - h) - (1 - h) * sqrt(1 - (1 - h)^2)) / π
                end
                @assert 0 <= w <= 0.5 "$w"
                _v = samplemesh(searchers, _inbound(origin - T(1.1) * d * n), vals, bg)
                w, _v, -n[1:N] / Z * sign(d)
            end
            filter!(!isnothing, v)

            edge_wts = getindex.(v, 1)
            edge_wt = sum(edge_wts)
            if edge_wt < 1
                wi = edge_wts / edge_wt
                center_normal = sum(wi .* getindex.(v, 3))
                Z = norm(center_normal)
                if Z > 0
                    center_normal /= Z
                    vedge = sum(wi .* getindex.(v, 2))
                end
            end
        end

        if isnothing(vedge)
            vcenter = samplemesh(searchers, origin, vals, bg)
            isPEC(vcenter) && return tensor ? 0 : vcenter
            return tensor ? (i == j) / vcenter : vcenter
        end

        length(center_normal) == 2 && (center_normal = [center_normal[1], center_normal[2], 0])
        origin = origin .- tol * center_normal
        vcenter = samplemesh(searchers, origin, vals, bg)
        isPEC(vcenter) && return tensor ? 0 : vcenter
        isPEC(vedge) && return tensor ? (i == j) / vcenter : vcenter
        # (isnothing(vedge) || (vcenter == vedge)) && return tensor ? (i == j) / vcenter : vcenter
        # ct2 += 1

        if tensor
            Pij = center_normal[i] * center_normal[j]
            # P = n * n'
            # isPEC(vedge) && return Pij * (1 - edge_wt) / vcenter
            # isPEC(vedge) && return Pij / vcenter
            Pij * ((1 - edge_wt) / vcenter + edge_wt / vedge) + ((i == j) - Pij) / ((1 - edge_wt) * vcenter + edge_wt * vedge)
            # P = improj(a)
            # P * mean(1 ./ a) + (In - P) / mean(a)
        else
            (1 - edge_wt) * vcenter + edge_wt * vedge
        end
    end .|> T
end

_inbound(v) = inbound(v, hcat(lb, ub), boundaries)

# if subpixel_smoothing
#     v = tmap(1:N) do i
#         tmap(1:i) do j
#             offset = (offsets("E$("xyz"[i])") + offsets("E$("xyz"[j])")) / 2
#             @time "subpixel smoothing entry" _supersamplemesh(T, N, searchers, vals, bg, offset, centers, Δs, _inbound, tensor, z, i, j)
#         end
#     end
#     return [j <= i ? v[i][j] : v[j][i] for i = 1:N, j = 1:N]
# end

# @unpack tree, tris, tri2body = searchers
# I, D = knn(tree, point, 4)
# tol = mean(D) / 100
# T = tris[I]
# S = ((point,) .- getindex.(T, 1)) .⋅ getindex.(T, 2)
# IS = zip(I, S)
# # IS = sort(collect(zip(I, S)), by=x -> abs(x[2]))
# io = Dict()
# for (i, s) in IS
#     b = tri2body[i]
#     if haskey(io, b)
#         io[b] = io[b] && (s < tol)
#     else
#         io[b] = s < tol
#     end
# end
# for b = keys(io)
#     io[b] && return vals[b]
# end
# bg


# function Searcher(meshes)
#     @time "sample trees" v = tmap(meshes) do m
#         tmap(collect(m)) do f
#             verts = reduce(hcat, vec3.(vertices(f)))
#             c = vec3(centroid(f))
#             n = vec3(ustrip(normal(f)))
#             (c, n, verts)
#         end |> filter(!isnothing)
#     end
#     tris = reduce(vcat, v)
#     tree = BallTree(reduce(hcat, getindex.(tris, 1)))
#     tri2body = reduce(vcat, fill.(eachindex(v), length.(v)))
#     Searcher(tree, tris, tri2body)
# end

function inbound(v, bbox, boundaries)
    tol = 1f-6
    N = size(bbox, 1)
    r = map(v[1:N], eachrow(bbox), eachrow(boundaries)) do v, (vl, vu), (bl, bu)
        if v < vl
            if bl == :PML || bl == :PEC
                return vl + tol
            elseif bl == :PMC || bl == :PEC
                return 2vl - v
            elseif bl == :periodic
                return vu + v - vl
            else
                error("Invalid boundary condition")
            end
        elseif v > vu
            if bu == :PML || bu == :PEC
                return vu - tol
            elseif bu == :PMC || bu == :PEC
                return 2vu - v
            elseif bu == :periodic
                return vl + v - vu
            else
                error("Invalid boundary condition")
            end
        else
            v
        end
    end
    length(v) > length(r) && return (r..., v[3])
    r
end

function raytri(ray, verts)
    try
        length(unique(sign.(verts \ ray))) == 1
    catch SingularException
        # @debug verts
        return raytri(ray, verts .+ ray / 100)
    end
end
