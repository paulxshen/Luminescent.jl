function refactor(d)
    N = size(d(1), 1)
    [namedtuple([k => d[k][i, :] for k = sort(keys(d))]) for i = 1:N]
end
_pmlfracs(a, N) = a[1:N, :]
_pmlfracs(a::Real, N) = fill(a, (N, 2))
_pmlfracs(a::AbstractVector, N) = stack([a[1:N], a[1:N]])

"""
    function setup(boundaries, sources, monitors, L, dx,                      approx_2D_mode=nothing; F=Float32)

Args
...
- L: vector of lengths in wavelengths of simulation domain
-                      approx_2D_mode: only applies to 2d which can be :TM (Ez, Hx, Hy) or :TE (Hz, Ex, Ey)
"""
function setup(dl, boundaries, sources, monitors, deltas, mode_deltas;
    approx_2D_mode=nothing,
    Ttrans=nothing, Tss=nothing, Tssmin=nothing,
    ϵ=1, μ=1, σ=0, m=0, γ=0, β=0,
    ϵ3=ϵ,
    F=Float32,
    Courant=0.9,
    deltas3=deltas,
    array=Array,
    pmlfracs=1,
    TEMP="",
    λ=1,
)

    if !isnothing(approx_2D_mode)
        approx_2D_mode = Symbol(approx_2D_mode)
    end
    N = length(deltas)
    pmlfracs = _pmlfracs(pmlfracs, N)
    (deltas, mode_deltas, ϵ, μ, σ, m, γ, β) = F.((deltas, mode_deltas, ϵ, μ, σ, m, γ, β))

    dx = minimum(minimum.(deltas))
    nres = int(1 / dx)

    mode_spacing = int(mode_deltas[1] / dl)
    L = size(ϵ) * dl
    sz0 = Tuple([isa(d, Number) ? int(l / d) : length(d) for (d, l) = zip(deltas, L)])
    sz = collect(sz0)

    # @show sz0,prod(sz0)

    _geometry = (; ϵ, μ, σ, m, γ, β) |> pairs |> OrderedDict
    geometry = OrderedDict()
    for (k, v) = pairs(_geometry)
        if isa(v, AbstractArray)
            # downsample(v, int(deltas / dl)) |> F
        elseif isa(v, Number)
            v = F(v)
            if k ∈ (:σ, :m)
                geometry[k] = fill(v, sz0)
            else
                geometry[k] = v
            end
            delete!(_geometry, k)
        end
    end
    _ϵ3 = ϵ3

    ϵmin, ϵmax = extrema(abs, _geometry.ϵ)
    μmin, μmax = extrema(abs, geometry.μ)
    nmax = sqrt(ϵmax * μmax)
    nmin = sqrt(ϵmin * μmin)

    dt = nmin / √(sum(v -> 1 / minimum(v)^2, deltas)) * Courant
    dt = 1 / ceil(1 / dt) |> F

    maxdeltas = maximum.(deltas)

    v = 0.15 / dt
    δ = -3log(1e-4) / nmin / 2 / (2v) |> F
    σpml = ϵmin * v
    mpml = μmin * v

    # r = max(1, maximum(16 ./ δ .* maxdeltas))
    # δ *= r
    # σpml /= r
    # m /= r

    pml_depths = max.(δ * pmlfracs, maxdeltas)
    pml_depths = ceil.(pml_depths ./ maxdeltas) .* maxdeltas
    @show σpml
    @show pml_depths

    if N == 1
        field_names = (:Ez, :Hy)
    elseif N == 2
        if approx_2D_mode == :TM
            Enames = (:Ez,)
            Hnames = (:Hx, :Hy)
            field_names = (:Ez, :Hx, :Hy)
        elseif approx_2D_mode == :TE
            Enames = (:Ex, :Ey)
            Hnames = (:Hz,)
            field_names = (:Ex, :Ey, :Hz)
        end
    else
        Enames = (:Ex, :Ey, :Ez)
        Hnames = (:Hx, :Hy, :Hz)
        field_names = (:Ex, :Ey, :Ez, :Hx, :Hy, :Hz)
    end

    nodes = fill(:U, N, 2)
    # v=Vector{Any}(fill(nothing, N))
    # field_boundvals = DefaultDict(() -> [copy(v),copy(v)])
    # geometry_padvals = DefaultDict(() -> [copy(v),copy(v)])
    db = Any[PML(j * i, pml_depths[i, int(j / 2 + 1.5)], σpml, mpml) for i = 1:N, j = (-1, 1)]
    field_boundvals = DefaultDict(() -> Array{Any,2}(fill(nothing, N, 2)))
    geometry_padvals = DefaultDict(() -> Array{Any,2}(fill(nothing, N, 2)))
    padamts = zeros(Int, N, 2)
    geometry_padamts = DefaultDict(() -> zeros(Int, N, 2))
    _geometry_padamts = DefaultDict(() -> zeros(Int, N, 2))
    is_field_on_lb = OrderedDict([k => zeros(Int, N) for k = field_names])
    is_field_on_ub = OrderedDict([k => zeros(Int, N) for k = field_names])
    bbox = zeros(F, N, 2)
    bbox[:, 2] .= L
    field_sizes = OrderedDict([k => collect(sz0) for k = field_names])

    for b = boundaries
        for i = b.dims
            if typeof(b) == Periodic
                db[i, :] = [Periodic(-abs(i)), Periodic(abs(i))]
            else
                if i > 0
                    db[i, 2] = b
                else
                    db[abs(i), 1] = b

                end
            end
        end
    end

    for i = 1:N
        for j = 1:2
            b = db[i, j]
            t = typeof(b)
            if t == PML
            elseif t == PEC
                nodes[i, j] = :E
            elseif t == PMC
                nodes[i, j] = :H
            end
        end
    end

    for i = 1:N
        if nodes[i, :] == [:U, :E]

            nodes[i, 1] = :H
        elseif nodes[i, :] == [:U, :H]

            nodes[i, 1] = :E
        elseif nodes[i, :] == [:E, :U]

            nodes[i, 2] = :H
        elseif nodes[i, :] == [:H, :U]

            nodes[i, 2] = :E
        elseif nodes[i, :] == [:U, :U]

            nodes[i, :] = [:E, :H]
        elseif nodes[i, :] == [:U, :U]

            nodes[i, :] = [:E, :H]
        elseif nodes[i, :] == [:E, :E]

        elseif nodes[i, :] == [:H, :H]

        end

    end


    for i = 1:N
        select = i .== 1:N
        xyz = para = perp = [:x, :y, :z]
        perp = [popat!(para, i)]
        for j = 1:2
            b = db[i, j]
            if isa(b, PML)
                npml = round(b.d / maxdeltas[i])
                _npml = npml * int(maxdeltas[i] / dl)
                # l1 = max.(1, round.(b.ramp_frac * l))
                # r1 = max.(1, round.(b.ramp_frac * r))
                # if any(l1 .> 0) || any(r1 .> 0)
                #     # rr=ReplicateRamp(convert.(F,b.σ))
                #     rr = convert.(F, b.σ)
                #     push!(geometry_padvals[:σ], OutPad(rr, l1, r1, sz0))
                #     push!(geometry_padvals[:m], OutPad(rr, l1, r1, sz0))
                # end
                if !isa(deltas[i], Number)
                    if j == 1
                        pushfirst!(deltas[i], fill(maxdeltas[i], npml)...)
                    else
                        push!(deltas[i], fill(maxdeltas[i], npml)...)
                    end
                end

                padamts[i, j] = npml
                for k = (:σ, :m, :ϵ, :μ)
                    geometry_padamts[k][i, j] = npml
                    _geometry_padamts[k][i, j] = _npml
                end
                for k = (:ϵ, :μ)
                    geometry_padvals[k][i, j] = :replicate
                end
                geometry_padvals[:σ][i, j] = TanhRamp(b.σ)
                geometry_padvals[:m][i, j] = TanhRamp(b.m)

                if j == 1
                    bbox[i, :] .-= b.d
                end

                sz[i] += npml
                for k = keys(field_sizes)
                    field_sizes[k][i] += npml
                end
            end
            f = nodes[i, j]
            for k = field_names
                q = startswith(String(k), String(f))
                if (q ? k[2] in para : k[2] in perp)
                    if isa(b, Periodic)
                        field_boundvals[k][i, j] = :periodic
                    else
                        field_boundvals[k][i, j] = 0
                    end
                end
            end
        end
    end


    for (k, v) = pairs(field_boundvals)
        is_field_on_lb[k] = !isnothing.(v[:, 1])
        is_field_on_ub[k] = !isnothing.(v[:, 2])
    end

    add_current_keys!(field_sizes)

    field_sizes = NamedTuple([k => Tuple(field_sizes[k]) for (k) = keys(field_sizes)])
    if N == 1
    elseif N == 3
        u0 = dict([k => zeros(F, Tuple(field_sizes[k])) for k = (:Ex, :Ey, :Ez, :Hx, :Hy, :Hz)])
    else
        if approx_2D_mode == :TM
            u0 = dict([k => zeros(F, Tuple(field_sizes[k])) for k = (:Ez, :Hx, :Hy)])
        else
            u0 = dict([k => zeros(F, Tuple(field_sizes[k])) for k = (:Ex, :Ey, :Hz)])
        end
    end
    add_current_keys!(u0)
    u0 = groupkeys(u0)

    geometry_sizes = NamedTuple([k => sz0 .+ sum(geometry_padamts[k], dims=2) for k = keys(geometry_padamts)])
    field_lims = OrderedDict{Symbol,Any}()
    for k = keys(geometry_sizes)
        if k in (:μ, :m)
            names = Hnames
        elseif k in (:ϵ, :σ)
            names = Enames
        end
        v = NamedTuple([
            begin
                xyz = f[2]
                terminations = zip(is_field_on_lb[f], is_field_on_ub[f])
                g = Symbol("$(k)$xyz$xyz")
                start = 0.5 - 0.5is_field_on_lb[f]
                v = [(start) (start + geometry_sizes[k] - 1)]
                f => v
            end for f = names
        ])
        field_lims[k] = v
    end
    field_lims = dict(merge(values(field_lims)...))
    add_current_keys!(field_lims)
    lb = bbox[:, 1]

    field_deltas = [_make_field_deltas(d, N, field_boundvals, field_sizes, i) for (i, d) = enumerate(deltas)]
    field_diffdeltas = [_make_field_deltas(d, N, field_boundvals, field_sizes, i, true) for (i, d) = enumerate(deltas)]

    geometry_padvals = NamedTuple(geometry_padvals)
    field_boundvals = NamedTuple(field_boundvals)
    field_diffpadvals = kmap(a -> reverse(a, dims=2), field_boundvals)
    field_spacings = int(field_deltas / dl)
    spacings = int(deltas / dl)
    field_diffpadvals = refactor(field_diffpadvals)

    sz = Tuple(sz)
    sz0 = Tuple(sz0)

    grid = (; F, N, L, bbox, sz0, sz, padamts, deltas, deltas3, lb, field_lims, field_sizes, field_boundvals, field_deltas, field_diffdeltas, field_diffpadvals, geometry_padvals, geometry_padamts, _geometry_padamts, dl, spacings, mode_spacing, mode_deltas,) |> dict

    println("making sources...")
    mode_solutions = []
    source_instances = SourceInstance.(sources, (grid,), (_ϵ3,), (TEMP,), (mode_solutions,))
    println("making monitors...")
    monitor_instances = MonitorInstance.(monitors, (grid,), (_ϵ3,), (TEMP,), (mode_solutions,))
    ϵeff = nothing
    # ϵeff = source_instances[1].ϵeff



    # geometry_padvals[:invϵ] = geometry_padvals[:ϵ]

    if N == 2
        # _geometry[:ϵ] = map(_geometry.ϵ) do x
        #     k, v = ϵeff
        #     k, v = F((k, v))
        #     abs(x - k) < 0.02 ? x - k + v : x
        # end
        geometry[:ϵ] = downsample(_geometry.ϵ, int(deltas / dl))
    end

    if !isa(Ttrans, Real)
        if endswith(string(Ttrans), "x")
            Ttrans = parse(F, Ttrans[1:end-1])
        else
            Ttrans = 1
        end
        Ttrans *= sum(L * nmax)
    end

    if Tss == nothing
        if isnothing(Tssmin)
            Tssmin = 120nmax / nres / N
        end
        v = reduce(vcat, wavelengths.(monitor_instances))
        v = Base.round.(v, digits=3)
        v = v |> Set |> collect |> sort |> reverse
        @show v
        if length(v) == 1
            Tss = 1
        else
            Tss = 1 / minimum(diff([0, (1 ./ v)...]))
            # Tss = ceil(100 / T) * T
        end
    end
    if !isnothing(Tssmin)
        @show Tssmin
        Tss *= Base.ceil(Int, Tssmin / Tss)
    end

    @show Ttrans, Tss
    Ttrans, Tss = convert.(F, (Ttrans, Tss))
    global prob = (;
                      grid,
                      source_instances, monitor_instances, field_names, approx_2D_mode, Courant,
                      Ttrans, Tss,
                      geometry, _geometry, nmax, nmin, ϵeff,
                      is_field_on_lb, is_field_on_ub,
                      u0, dt, array, λ) |> pairs |> OrderedDict

    _gpu = x -> gpu(array, x)
    if array == Array
        println("using CPU backend.")
    else
        println("using GPU backend.")
    end
    for k = keys(prob)
        if k in (:u0, :_geometry, :geometry, :source_instances, :monitor_instances)
            prob[k] = _gpu(prob[k],)
        end
    end
    prob.grid[:field_diffdeltas] = _gpu.(prob.grid[:field_diffdeltas])
    prob
end
update = update
setup = setup