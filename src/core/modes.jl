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
    @nograd deltas
    p = 0
    p += conj(u(:Ex, 0)) ⊙ v(:Hy, 0)
    p -= conj(u(:Ey, 0)) ⊙ v(:Hx, 0)
    p += conj(u(:Hy, 0)) ⊙ v(:Ex, 0)
    p -= conj(u(:Hx, 0)) ⊙ v(:Ey, 0)
    sum(reduce(deltas; init=p) do a, v
        a .* v
    end)
end

function normalize_mode(m, deltas)
    p = inner(m, m, deltas)
    p = real(p)
    m = m / sqrt(abs(p))
    if p < 0
        mirror_mode(m)
    else
        m
    end
end
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
    m = kmap(m) do a
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


function solvemodes(ϵ, dl, λ, neigs, spacing, path; mode_solutions=nothing)
    if !isnothing(mode_solutions)
        i = findfirst(mode_solutions) do v
            _ϵ, _dl, _λ, _neigs = v
            size(ϵ) == size(_ϵ) && sum(abs, ϵ - _ϵ) / sum(abs, ϵ) < 1e-3 && dl == _dl && λ == _λ && neigs <= _neigs
        end
        if !isnothing(i)
            return mode_solutions[i][end][1:neigs]
        end
    end

    # m = round((size(ϵ, 1) - size(ϵ, 2)) / 2)
    # n = size(ϵ, 1) - size(ϵ, 2) - m
    # ϵ = pad(ϵ, :replicate, [0, m], [0, n])
    name = round(1000000Float64(λ))

    if !isfile(joinpath(path, "$(name)_mode_$(neigs-1).npz"))
        println("run empy")
        npzwrite(joinpath(path, "args.npz"), Dict("eps" => begin
                ϵ1 = (ϵ[1:end-1, :] + ϵ[2:end, :]) / 2
                (ϵ1[:, 1:end-1] + ϵ1[:, 2:end]) / 2
            end,
            "dl" => dl, "center_wavelength" => λ, "neigs" => neigs, "name" => name))
        fn = joinpath(path, "solvemodes.py")
        try
            run(`python $fn $path`)
        catch e
            run(`python3 $fn $path`)
        end
    end

    modes = [npzread(joinpath(path, "$(name)_mode_$(i-1).npz")) for i = 1:neigs]
    modes = [merge(mode, OrderedDict(["J$s" => mode["E$s"] .* ϵ for s = "xy"])) for mode in modes]
    global modes = [SortedDict([Symbol(k) => downsample(mode(k), spacing) for k = keys(mode) if string(k)[end] in "xy"]) |> pairs |> NamedTuple for mode in modes]

    if !isnothing(mode_solutions)
        println("saving mode solutions")
        push!(mode_solutions, (ϵ, dl, λ, neigs, modes))
    end
    modes
end