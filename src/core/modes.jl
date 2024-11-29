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

# function normalize_mode(m, deltas)
function inner(u, v, deltas)
    p = 0
    if haskey(u, :Ex) && haskey(v, :Hy)
        p += conj(u.Ex) .* v.Hy
    end
    if haskey(u, :Ey) && haskey(v, :Hx)
        p -= conj(u.Ey) .* v.Hx
    end
    if haskey(v, :Ex) && haskey(u, :Hy)
        p += conj(u.Hy) .* v.Ex
    end
    if haskey(v, :Ey) && haskey(u, :Hx)
        p -= conj(u.Hx) .* v.Ey
    end
    sum(p .* ignore_derivatives() do
        reduce(.*, deltas) / 2
    end)
end


function keepxy(mode)
    namedtuple([k => mode[k] for k in keys(mode) if any(endswith.((string(k),), ("x", "y")))])
end


function mode_decomp(m, u, deltas)
    p = real(inner(m, m, deltas))
    m1 = m / ignore_derivatives() do
        sqrt(abs(p))
    end
    m2 = mirror_mode(m1)
    if p < 0
        m2, m1 = m1, m2
    end
    # @ignore_derivatives_vars m1, m2, deltas
    inner.((m1, m2), (u,), (deltas,))
end


function collapse_mode(m, Eonly=false, p=:TE)
    # @unpack Ex, Ey, Ez, Hx, Hy, Hz = m
    # Ex, Ey, Ez, Hx, Hy, Hz = map([Ex, Ey, Ez, Hx, Hy, Hz]) do a
    #     mean(a, dims=2) |> vec
    # end
    m = kmap(m) do a
        mean(a, dims=2) |> vec
    end
    # (; Ex, Hy, Ez), sum(E .* ϵ, dims=2) ./ sum(E, dims=2) |> vec
    if p == :TE
        Eonly ? (; Ex=m.Ex,) : (; Ex=m.Ex, Hy=m.Hy)
        #, maximum.(eachrow(ϵ))
    else
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

Base.convert(::Type{Float64}, x::ComplexF64) = real(x)

function solvemodes(ϵ, dl, λ, neigs, spacing, path)
    # m = round((size(ϵ, 1) - size(ϵ, 2)) / 2)
    # n = size(ϵ, 1) - size(ϵ, 2) - m
    # ϵ = pad(ϵ, :replicate, [0, m], [0, n])
    @show λ, dl
    npzwrite(joinpath(path, "args.npz"), Dict("eps" => imresize(ϵ, size(ϵ) - 1), "dl" => dl, "wl" => λ, "neigs" => neigs))
    fn = joinpath(path, "solvemodes.py")
    run(`python $fn $path`)
    modes = [npzread(joinpath(path, "mode$(i-1).npz")) for i = 1:neigs]
    # modes = [NamedTuple([Symbol(k) => v[:, m+1:end-n] for (k, v) in mode]) for mode in modes]
    display(heatmap(ϵ))
    display(heatmap(real(modes[1].Ex)))
    display(heatmap(real(modes[1].Hy)))
    modes = [NamedTuple([Symbol(k) => downsample(mode(k), spacing) for k = (:Ex, :Ey, :Hx, :Hy)]) for mode in modes]
    modes
end
#     x = range(dl / 2; step=dl, length=size(ϵ, 1))
#     y = range(dl / 2; step=dl, length=size(ϵ, 2))
#     # x = range(0; step=dl, length=size(ϵ, 1) + 1)
#     # y = range(0; step=dl, length=size(ϵ, 2) + 1)
#     tol = 1e-8
#     boundary = (0, 0, 0, 0)
#     f = (x, y) -> begin
#         v = getindexf(ϵ, x / dl + 0.5, y / dl + 0.5)
#         (v, 0, 0, v, v)
#     end
#     # global _as = ϵ, x, y
#     solver = VectorModesolver.VectorialModesolver(λ, x, y, boundary, f)
#     modes = VectorModesolver.solve(solver, neigs, tol)
#     # plot_mode_fields(modes[1]) |> display
#     # error()
#     # display(heatmap(ϵ))
#     # display(heatmap(real(transpose(modes[1].Ex))))
#     # for i = 1:2
#     #     Ex = modes[i].Ex
#     #     e = mean(Ex .* ϵ, dims=2) ./ mean(Ex, dims=2)
#     #     display(lines(abs.(vec(e))))
#     # end
#     # error()
#     modes = [namedtuple([k => getfield(mode, k)[:, m+1:end-n] for k = (:Ex, :Ey, :Hx, :Hy)]) for mode in modes]


#     # [namedtuple([k => transpose(getfield(mode, k)) for k = (:Ex, :Ey, :Hx, :Hy)]) for mode in modes]
# end