function normalize_mode(m, dx)
    @unpack Ex, Hy = m
    @show p = real(Ex ⋅ Hy) * dx^ndims(Ex) / 2
    # (; Ex=Ex / √p, Hy=Hy / √p, Ez=Ez / √p)
    (; Ex=Ex / √p, Hy=Hy / √p)
end

function mode_decomp(m, u,)
    Ex = m.Ex ⋅ field(u, :Ex) / norm(m.Ex)^2
    Hy = m.Hy ⋅ field(u, :Hy) / norm(m.Hy)^2
    ap = Hy + Ex
    am = Hy - Ex
    [ap, am]
end


function collapse_mode(m, ϵ)
    @unpack Ex, Ey, Ez, Hx, Hy, Hz = m
    E = sqrt.(abs2.(Ex) + abs2.(Ey) + abs2.(Ez))
    Ex, Ez, Hy = map([Ex, Ez, Hy]) do a
        mean(a, dims=2) |> vec
    end
    (; Ex, Hy, Ez), sum(E .* ϵ, dims=2) ./ sum(E, dims=2) |> vec
end

function reframe(frame, u, inv=false)
    d = p = signs = sz = invdims = 0
    ignore_derivatives() do
        d = ndims(u[1]) + 1
        p = findfirst.(x -> abs(1 - abs(x)) < 1.0f-3, frame)
        signs = sign.(getindex.(frame, p))
        invdims = findall(isequal(-1), signs) |> Tuple
        # println(p)
    end
    if inv
        ignore_derivatives() do
            sz = collect(size(u[1]))
            insert!(sz, p[3], 1)
            if d == 2
                insert!(sz, p[2], 1)
                (@assert p[2] == 3)
            end
            sz = Tuple(sz)
        end
        u = [reshape(a, sz) for a in u]
        u = permutedims.(u, (p,))
        u = u[p]

        u = [s * reverse(a, dims=invdims) for (s, a) in zip(signs, u)]
        u = dropdims.(u, dims=Tuple(d:3))
        global q = u
        # p = invperm(p)
    else
        u = [reshape(a, size(a)..., fill(1, 3 - ndims(a))...) for a in u]
        # invpermute!(u, perm)
        u = [s * reverse(a, dims=invdims) for (s, a) in zip(signs, u)]
        u = [u[i] for i = invperm(p)]
        u = permutedims.(u, (invperm(p),))
    end
end
invreframe(frame, u) = reframe(frame, u, true)