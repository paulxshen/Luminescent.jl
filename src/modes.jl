mirror_mode(m) = namedtuple([k =>
    begin
        a = reverse(m[k], dims=1)
        v = ignore_derivatives() do
            string(k)[2:2]
        end
        if v in ("x", "z")
            a = -a
        end
        a
    end for k in keys(m)])

# function normalize_mode(m, dx)
function inner(u, v; dx=1)
    p = 0
    if haskey(u, :Ex) && haskey(v, :Hy)
        p += u.Ex ⋅ v.Hy
    end
    if haskey(u, :Ey) && haskey(v, :Hx)
        p -= u.Ey ⋅ v.Hx
    end
    if haskey(v, :Ex) && haskey(u, :Hy)
        p += u.Hy ⋅ v.Ex
    end
    if haskey(v, :Ey) && haskey(u, :Hx)
        p -= u.Hx ⋅ v.Ey
    end
    p * ignore_derivatives() do
        dx^ndims(first(values(u))) / 2
    end
    # if p > 0
    #     return r
    # else
    #     return mirror_mode(r)
    # end
end


function keepxy(mode)
    namedtuple([k => mode[k] for k in keys(mode) if any(endswith.((string(k),), ("x", "y")))])
end


function mode_decomp(m, u, dx)
    p = real(inner(m, m; dx))
    m1 = m / sqrt(abs(p))
    m2 = mirror_mode(m1)
    if p < 0
        m2, m1 = m1, m2
    end
    # k = ignore() do
    #     intersect(keys(m1), keys(u)) |> collect |> sort
    # end
    # m = [vcat(getindex.((m,), k)...) for m = (m1, m2)]
    # u = vcat(getindex.((u,), k)...)
    inner.((m1, m2), (u,); dx)
end


function collapse_mode(m, p=:TE)
    @unpack Ex, Ey, Ez, Hx, Hy, Hz = m
    Ex, Ey, Ez, Hx, Hy, Hz = map([Ex, Ey, Ez, Hx, Hy, Hz]) do a
        mean(a, dims=2) |> vec
    end
    # (; Ex, Hy, Ez), sum(E .* ϵ, dims=2) ./ sum(E, dims=2) |> vec
    if p == :TE
        (; Ex, Hy, Ez)
        #, maximum.(eachrow(ϵ))
    else
        (; Hx, Ey, Hz)
    end
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