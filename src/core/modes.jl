

# mirror_mode(m) =
#     kvmap(m) do k, a
#         k = string(k)
#         if k[2:2] in ("x", "z")
#             a = reverse(a, dims=1)
#             if k[1] == 'E'
#                 a = -a
#             end
#         else
#             if k[1] == 'H'
#                 a = -a
#             end
#         end
#         a
#     end
# mirror_mode(m) =
#     kvmap(m) do k, a
#         a = reverse(a, dims=1)
#         if string(k)[2:2] in ("x", "z")
#             a = -a
#         end
#         k => a
#     end
mirror_mode(m) =
    kvmap(m) do k, a
        if string(k)[1] == 'H'
            a = -a
        end
        k => a
    end

function _inner(f, u, v)
    f.(sum(((:Ex, :Hy, 1), (:Ey, :Hx, -1), (:Hy, :Ex, 1), (:Hx, :Ey, -1))) do (i, j, s)
        conj(u(i, 0)) .* v(j, 0) * s
    end / 2)
end

inner(f::Function, u, v, dA) = sum(dA .* _inner(f, u, v))
inner(u, v, dA) = inner(identity, u, v, dA)
power(u, dA) = real(inner(u, u, dA))
function normalize_mode(m, dA)
    p = power(m, dA)
    if p < 0
        m = mirror_mode(m)
    end
    m / sqrt(abs(p))
end

poynting(u) = [u.Ey .* u.Hz - u.Ez .* u.Hy,
    u.Ez .* u.Hx - u.Ex .* u.Hz,
    u.Ex .* u.Hy - u.Ey .* u.Hx]
# function keepxy(mode)
#     namedtuple([k => mode[k] for k in keys(mode) if any(endswith.((string(k),), ("x", "y")))])
# end


# function mode_decomp(m, u, deltas)
#     @nograd deltas, m
#     p = inner(m, m, deltas)
#     @assert imag(p) < 1e-3
#     p = real(p)
#     m1 = m / sqrt(abs(p))
#     m2 = mirror_mode(m1)
#     if p < 0
#         m2, m1 = m1, m2
#     end
#     inner.((m1, m2), (u,), (deltas,))
# end


function collapse_mode(m, p=:TE)
    p = Symbol(p)
    m = vmap(m) do a
        mean(a, dims=2) |> vec
    end
    if p == :TE
        (; Ex=m.Ex, Hy=m.Hy, Jx=m.Jx)
    elseif p == :TM
        (; Ey=m.Ey, Hx=m.Hx, Jy=m.Jy)
    else
        error("invalid polarization")
    end
end

function insert(a, i, v)
    if i == 1
        [v, a...]
    else
        [a[1:i-1]..., v, a[i:end]...]
    end
end

function localframe(u, monitor::PlaneMonitorInstance; approx_2D_mode=nothing)
    @nograd monitor
    R = @ignore_derivatives inv(monitor.frame)
    Is = @ignore_derivatives monitor.local_plane_Is

    if isnothing(approx_2D_mode)
        u = todvoa(u)
        u = vmap(u) do v
            R * v
        end
        u = todxyz(u)
    else
        u = (E=R * [u.Ex, u.Ey], H=[u.Hz])
        u = (Ex=u.E[1], Hy=u.H[1])
        # Ex = getindexf.((u.Ex,), Is)
        # Hy = getindexf.((u.Hy,), Is)
        # Ez = getindexf.((u.Ez,), Is)
        # (; Ex, Hy, Ez)
    end
    a = vmap(u) do a
        a = map(Is) do I
            getindexf(a, I..., :)
        end
        n = length(a[end])
        a = reduce(vcat, a)
        a = reshape(a, n, size(Is)...)
        a = permutedims(a, ((2:ndims(a))..., 1))
    end
end
function localframe(u, m::SphereMonitorInstance; kwargs...)
    @unpack As = m
    nf = size(u(1))[end]
    np = length(m)
    u = kvmap(u) do k, v
        # k => As(k) * v
        k => stack(map(eachcol(v)) do v
            As(k) * v
        end)#temp
        # k => Array(As(k)) * Array(v) |> typeof(v) #temp
    end

    u = vmap(todvoa(u)) do u
        a = m.invframes
        # map(, zip(u...)) do R, v
        #     R * v
        # end
        # map(axes(a, ndims(a)), eachrow.(u)...) do i, v...
        v = map(eachslice(a; dims=ndims(a)), eachrow.(u)...) do a, v...
            a * stack(v)'
        end
        permutedims(stack(v), (3, 2, 1))
    end
    u = todxyz(vmap(u) do v
        reshape.(eachslice(v; dims=ndims(v)), ((size(m)..., nf),))
    end)
end


addJ(mode, ϵ) = NamedTuple([Symbol("J$s") => mode("E$s", 0) .* ϵ .* .!isPEC(ϵ) for s = "xyz"])
