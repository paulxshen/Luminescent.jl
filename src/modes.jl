function normalize_mode(m, dx)
    @unpack Ex, Hy, = m
    # @unpack Ex, Hy, Ez = m
    p = real(Ex ⋅ Hy) * dx^ndims(Ex) / 2
    m / √p
    # (; Ex=Ex / √p, Hy=Hy / √p, Ez=Ez / √p)
end

function keepxy(mode)
    dict([k => mode[k] for k in keys(mode) if string(k)[2] ∈ ('x', 'y')])
end


# function mode_decomp(m, u::AbstractVector, dx)
#     Ex,Hy,Hz
function mode_decomp(m, u, dx)

    polarization = get_polarization(u)
    m = normalize_mode(m, dx)
    if polarization == :TM
        ap_TE = am_TE = 0
    else
        Ex = m.Ex ⋅ field(u, :Ex) / norm(m.Ex)^2
        Hy = m.Hy ⋅ field(u, :Hy) / norm(m.Hy)^2
        ap_TE = Hy + Ex
        am_TE = Hy - Ex
    end
    if polarization == :TE
        ap_TM = am_TM = 0
    else
        Hx = m.Hx ⋅ field(u, :Hx) / norm(m.Hx)^2
        Ey = m.Ey ⋅ field(u, :Ey) / norm(m.Ey)^2
        am_TM = Ey + Hx
        ap_TM = Ey - Hx
    end
    # [ap_TE + ap_TM, am_TE + am_TM]
    if polarization == :TM
        return [ap_TM, am_TM]
    end
    if polarization == :TE
        return [ap_TE, am_TE]
    end
    return [ap_TE, am_TE]
end


function collapse_mode(m, ϵ)
    @unpack Ex, Ey, Ez, Hx, Hy, Hz = m
    E = sqrt.(abs2.(Ex) + abs2.(Ey) + abs2.(Ez))
    Ex, Ez, Hy = map([Ex, Ez, Hy]) do a
        mean(a, dims=2) |> vec
    end
    # (; Ex, Hy, Ez), sum(E .* ϵ, dims=2) ./ sum(E, dims=2) |> vec
    (; Ex, Hy, Ez), maximum.(eachrow(ϵ))
end

function insert(a, i, v)
    if i == 1
        [v, a...]
    else
        [a[1:i-1]..., v, a[i:end]...]
    end
end
function reframe(frame, u, inv=false)
    # p = = 0
    d = ndims(u[1]) + 1
    p = findfirst.(x -> abs(1 - abs(x)) < 1.0f-3, frame)
    signs = sign.(getindex.(frame, p))
    invdims = findall(isequal(-1), signs) |> Tuple
    if inv
        sz = size(u[1])
        sz = insert(sz, p[3], 1)
        if d == 2
            sz = insert(sz, p[2], 1)
            # (@assert p[2] == 3)
        end
        sz = Tuple(sz)

        u = [reshape(a, sz) for a in u]
        u = permutedims.(u, (p,))
        u = u[p]

        u = [s * reverse(a, dims=invdims) for (s, a) in zip(signs, u)]
        u = dropdims.(u, dims=Tuple(d:3))
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