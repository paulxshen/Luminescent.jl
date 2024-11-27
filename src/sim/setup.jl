"""
    function setup(boundaries, sources, monitors, L, dx, polarization=nothing; F=Float32)

Args
...
- L: vector of lengths in wavelengths of simulation domain
- polarization: only applies to 2d which can be :TM (Ez, Hx, Hy) or :TE (Hz, Ex, Ey)
"""
function setup(dl, boundaries, sources, monitors, deltas, mode_deltas;
    polarization=:TE,
    transient_duration=0, steady_state_duration=0,
    ϵ=1, μ=1, σ=0, m=0,
    ϵ3=1,
    F=Float32,
    pml_depths=nothing, pml_ramp_fracs=0.2,
    Courant=0.9,
    deltas3=deltas,
    kw...)
    (deltas, mode_deltas, ϵ, μ, σ, m) = F.((deltas, mode_deltas, ϵ, μ, σ, m))
    N = length(deltas)
    L = size(ϵ) * dl
    sz = Tuple([isa(d, Number) ? int(l / d) : length(d) for (d, l) = zip(deltas, L)])
    a = ones(F, Tuple(sz))
    _geometry = (; ϵ, μ, σ, m) |> pairs |> OrderedDict
    geometry = OrderedDict()
    # global _a = geometry, _geometry, ϵ3
    for (k, v) = pairs(_geometry)
        geometry[k] = if isa(v, AbstractArray)
            @time downsample(v, int(deltas / dl))
        elseif k ∈ (:σ, :m)
            a * v
        else
            v
        end
    end
    ϵ3 = if N == 3
        geometry.ϵ
    else
        downsample(ϵ3, int(deltas3 / dl))
    end

    ϵmin, ϵmax = extrema(geometry.ϵ)
    μmin, μmax = extrema(geometry.μ)
    nmax = sqrt(ϵmax * μmax)
    nmin = sqrt(ϵmin * μmin)
    N = length(sz)
    if isa(pml_depths, Number)
        pml_depths = fill(pml_depths, N)
    end
    # if isa(pml_ramp_fracs, Number)
    #     pml_ramp_fracs = fill(pml_ramp_fracs, N)
    # end
    maxdeltas = maximum.(deltas)
    if isnothing(pml_depths)
        mpml = σpml = 2.0
        @show δ = 4 / nmin / (σpml + mpml)
        pml_depths = max.([δ, δ, 0.2δ][1:N], maxdeltas)
    end
    pml_depths = trim.(pml_depths, maxdeltas)



    if N == 1
        field_names = (:Ez, :Hy)
        polarization = nothing
    elseif N == 2
        if polarization == :TM
            Enames = (:Ez,)
            Hnames = (:Hx, :Hy)
            field_names = (:Ez, :Hx, :Hy)
        elseif polarization == :TE
            Enames = (:Ex, :Ey)
            Hnames = (:Hz,)
            field_names = (:Ex, :Ey, :Hz)
        end
    else
        polarization = nothing
        Enames = (:Ex, :Ey, :Ez)
        Hnames = (:Hx, :Hy, :Hz)
        field_names = (:Ex, :Ey, :Ez, :Hx, :Hy, :Hz)
    end

    nodes = fill(:U, N, 2)
    # pml_depths = [xpml, ypml, zpml]
    # pml_ramp_fracs = [xpml_ramp_frac, ypml_ramp_frac, zpml_ramp_frac]
    db = Any[PML(j * i, pml_depths[i], σpml, mpml) for i = 1:N, j = (-1, 1)]
    field_boundvals = DefaultDict(() -> Array{Any,2}(fill(nothing, N, 2)))
    geometry_padvals = DefaultDict(() -> Array{Any,2}(fill(nothing, N, 2)))
    geometry_padamts = DefaultDict(() -> zeros(Int, N, 2))
    _geometry_padamts = DefaultDict(() -> zeros(Int, N, 2))
    is_field_on_lb = Dict([k => zeros(Int, N) for k = field_names])
    is_field_on_ub = Dict([k => zeros(Int, N) for k = field_names])
    bbox = zeros(F, N, 2)
    bbox[:, 2] .= L
    field_sizes = Dict([k => collect(sz) for k = field_names])

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
                _npml = round(b.d / dl)
                # l1 = max.(1, round.(b.ramp_frac * l))
                # r1 = max.(1, round.(b.ramp_frac * r))
                # if any(l1 .> 0) || any(r1 .> 0)
                #     # rr=ReplicateRamp(convert.(F,b.σ))
                #     rr = convert.(F, b.σ)
                #     push!(geometry_padvals[:σ], OutPad(rr, l1, r1, sz))
                #     push!(geometry_padvals[:m], OutPad(rr, l1, r1, sz))
                # end
                if !isa(deltas[i], Number)
                    if j == 1
                        pushfirst!(deltas[i], fill(maxdeltas[i], npml)...)
                    else
                        push!(deltas[i], fill(maxdeltas[i], npml)...)
                    end
                end
                for k = (:σ, :m, :ϵ, :μ)
                    geometry_padamts[k][i, j] = npml
                    _geometry_padamts[k][i, j] = _npml
                end
                for k = (:ϵ, :μ)
                    geometry_padvals[k][i, j] = :replicate
                end
                geometry_padvals[:σ][i, j] = b.σ
                geometry_padvals[:m][i, j] = b.m

                if j == 1
                    bbox[i, :] .-= b.d
                end
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
    fielmaxdeltas = NamedTuple([k => zeros(F, Tuple(field_sizes[k])) for (k) = field_names])
    if N == 1
    elseif N == 3
        u0 = NamedTuple([k => zeros(F, Tuple(field_sizes[k])) for k = (:Ex, :Ey, :Ez, :Hx, :Hy, :Hz, :Jx, :Jy, :Jz)])
    else
        if polarization == :TM
            u0 = NamedTuple([k => zeros(F, Tuple(field_sizes[k])) for k = (:Ez, :Hx, :Hy, :Jz)])
        else
            u0 = NamedTuple([k => zeros(F, Tuple(field_sizes[k])) for k = (:Ex, :Ey, :Hz)])
        end
    end

    geometry_sizes = NamedTuple([k => sz .+ sum(geometry_padamts[k], dims=2) for k = keys(geometry_padamts)])
    # field_lims = NamedTuple(Pair.(keys(geometry_sizes), Base.oneto.(values(geometry_sizes))))
    field_lims = OrderedDict{Symbol,Any}()
    field_grimaxdeltas = Dict{Symbol,Any}()
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
    field_lims = merge(values(field_lims)...)
    field_lims = add_current_keys(field_lims)
    lb = bbox[:, 1]

    field_deltas = [_make_field_deltas(d, N, field_boundvals, field_sizes, i) for (i, d) = enumerate(deltas)]
    field_diffdeltas = [_make_field_deltas(d, N, field_boundvals, field_sizes, i, true) for (i, d) = enumerate(deltas)]

    geometry_padvals = NamedTuple(geometry_padvals)
    field_boundvals = NamedTuple(field_boundvals)
    field_diffpadvals = kmap(a -> reverse(a, dims=2), field_boundvals)
    field_spacings = int(field_deltas / dl)
    spacings = int(deltas / dl)

    grid = (; F, N, L, bbox, sz, deltas, deltas3, lb, field_lims, field_sizes, field_boundvals, field_deltas, field_diffdeltas, field_diffpadvals, geometry_padvals, geometry_padamts, _geometry_padamts, dl, spacings,)

    source_instances = SourceInstance.(sources, (grid,), (ϵ3,))
    monitor_instances = MonitorInstance.(monitors, (grid,), (ϵ3,))

    dt = nmin / √(sum(dx -> 1 / minimum(dx)^2, deltas)) * Courant
    dt = 1 / ceil(1 / dt) |> F
    sz = Tuple(sz)

    # geometry_padvals[:invϵ] = geometry_padvals[:ϵ]



    if transient_duration == 0
        transient_duration = sum(L * nmax)
    end
    if steady_state_duration == 0
        if isempty(monitors)
            steady_state_duration = 1
        else
            v = reduce(vcat, wavelengths.(monitor_instances))
            v = v |> Set |> collect |> sort |> reverse
            if length(v) == 1
                steady_state_duration = 6
            else
                steady_state_duration = 6 / minimum(diff([0, (1 ./ v)...]))
            end
        end
    end
    transient_duration, steady_state_duration = convert.(F, (transient_duration, steady_state_duration))
    global res = (;
                     grid,
                     source_instances, monitor_instances, field_names,
                     mode_deltas,
                     polarization, Courant,
                     transient_duration, steady_state_duration,
                     geometry, _geometry, nmax, nmin,
                     is_field_on_lb, is_field_on_ub,
                     u0, dt, kw...) |> pairs |> OrderedDict

end
update = update
setup = setup