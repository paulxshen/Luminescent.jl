using Porcupine: keys, values
function get_polarization(u)
    # u=flatten(u)
    if haskey(u, :E)
        if length(u.E) == 2
            return :TE
        elseif length(u.H) == 2
            return :TM
        end
    else
        if haskey(u, :Hy)
            return :TE
        elseif haskey(u, :Ey)
            return :TM
        end
    end
    nothing
end
function add_current_keys(N)
    N = OrderedDict(pairs(N))
    add_current_keys!(N)
end
function add_current_keys!(N::AbstractDict)
    for k in keys(N) |> collect
        if startswith(String(k), "E")
            N[Symbol("J" * String(k)[2:end])] = N[k]
        end

        # if startswith(String(k), "H")
        #     N[Symbol("M" * String(k)[2:end])] = N[k]
        # end
    end
    N
end
"""
    function setup(boundaries, sources, monitors, L, dx, polarization=nothing; F=Float32)

Args
...
- L: vector of lengths in wavelengths of simulation domain
- polarization: only applies to 2d which can be :TM (Ez, Hx, Hy) or :TE (Hz, Ex, Ey)
"""
function setup(boundaries, sources, monitors, Δ, dl;
    polarization=:TE,
    # transient_duration=max_source_dist(sources), steady_state_duration=2,
    transient_duration=0, steady_state_duration=0,
    ϵ=1, μ=1, σ=0, m=0,
    F=Float32,
    pml_depths=nothing, pml_ramp_fracs=0.2,
    # xpml=0.4, ypml=0.4, zpml=0.4,
    # xpml_ramp_frac=0.5, ypml_ramp_frac=0.5, zpml_ramp_frac=0.5,

    # ϵmin=1,
    # Courant=0.5,
    Courant=0.95,
    kw...)
    ratio = Int.(Δ / dl)
    Δ, ϵ, μ, σ, m = [convert.(F, a) for a in (Δ, ϵ, μ, σ, m)]
    N = length(Δ)
    if isa(Δ, Number)
        Δ = fill(Δ, N)
    end
    L = size(ϵ) * dl
    sz = Tuple(round.(Int, L ./ Δ))
    a = ones(F, Tuple(sz))
    global _geometry = (; ϵ, μ, σ, m) |> pairs |> OrderedDict
    global geometry = OrderedDict()
    # f(x::Number) = a * x
    # f(x::AbstractArray) = x
    # geometry = kmap(f, geometry)
    for (k, v) = pairs(_geometry)
        if isa(v, AbstractArray)
            geometry[k] = imresize(v, sz)
        else
            geometry[k] = a * v
        end
    end

    ϵmin, ϵmax = extrema(geometry.ϵ)
    μmin, μmax = extrema(geometry.μ)
    @show nmax = sqrt(ϵmax * μmax)
    @show nmin = sqrt(ϵmin * μmin)
    N = length(sz)
    if isa(pml_depths, Number)
        pml_depths = fill(pml_depths, N)
    end
    # if isa(pml_ramp_fracs, Number)
    #     pml_ramp_fracs = fill(pml_ramp_fracs, N)
    # end
    if isnothing(pml_depths)
        mpml = σpml = 2.0
        @show δ = 4 / nmin / (σpml + mpml)
        pml_depths = max.([δ, δ, 0.2δ][1:N], Δ)
    end
    pml_depths = trim.(pml_depths, Δ)



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
                nspml = Int.(b.d / Δ[i])
                # l1 = max.(1, round.(b.ramp_frac * l))
                # r1 = max.(1, round.(b.ramp_frac * r))
                # if any(l1 .> 0) || any(r1 .> 0)
                #     # rr=ReplicateRamp(convert.(F,b.σ))
                #     rr = convert.(F, b.σ)
                #     push!(geometry_padvals[:σ], OutPad(rr, l1, r1, sz))
                #     push!(geometry_padvals[:m], OutPad(rr, l1, r1, sz))
                # end
                # l2 = l - l1
                # r2 = r - r1
                # if any(l2 .> 0) || any(r2 .> 0)

                for k = (:σ, :m, :ϵ, :μ)
                    geometry_padamts[k][i, j] = nspml
                end
                for k = (:ϵ, :μ)
                    geometry_padvals[k][i, j] = :replicate
                end
                geometry_padvals[:σ][i, j] = b.σ
                geometry_padvals[:m][i, j] = b.m

                if j == 1
                    bbox[i, :] .+= b.d
                end
                for k = keys(field_sizes)
                    field_sizes[k][i] += nspml
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
    fields = NamedTuple([k => zeros(F, Tuple(field_sizes[k])) for (k) = field_names])
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
    # fieldlims = NamedTuple(Pair.(keys(geometry_sizes), Base.oneto.(values(geometry_sizes))))
    global fieldlims = OrderedDict{Symbol,Any}()
    field_grids = Dict{Symbol,Any}()
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
        fieldlims[k] = v
    end
    fieldlims = merge(values(fieldlims)...)
    fieldlims = add_current_keys(fieldlims)
    origin = -bbox[:, 1]

    global _origin, _fieldlims, _sources = origin, fieldlims, sources
    source_instances = SourceInstance.(sources, (Δ,), (field_sizes,), (origin,), (fieldlims,); F)
    monitor_instances = MonitorInstance.(monitors, (Δ,), (origin,), (fieldlims,); F)
    # roi = MonitorInstance(Monitor(zeros(N), zeros(N), Δ * sz), Δ, sz, common_left_pad_amount, is_field_on_lb, fl; F)
    onedge = NamedTuple([k => hcat(v, is_field_on_ub[k]) for (k, v) = pairs(is_field_on_lb)])

    dt = nmin / √(sum(dx -> 1 / dx^2, Δ)) * Courant
    dt = 1 / ceil(1 / dt) |> F
    sz = Tuple(sz)

    # geometry_padvals[:invϵ] = geometry_padvals[:ϵ]
    geometry_padvals = NamedTuple(geometry_padvals)
    field_boundvals = NamedTuple(field_boundvals)

    if transient_duration == 0
        transient_duration = sum(sz .* Δ * nmax)
    end
    if steady_state_duration == 0
        if isempty(monitors)
            steady_state_duration = 1
        else
            v = reduce(vcat, wavelengths.(monitors))
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
                     geometry_padvals, geometry_padamts, field_boundvals, fieldlims, bbox,
                     source_instances, monitor_instances, field_names,
                     polarization, F, Courant,
                     transient_duration, steady_state_duration,
                     geometry, _geometry, nmax, nmin,
                     is_field_on_lb, is_field_on_ub, geometry_sizes,          # roi,
                     u0, fields, N, sz, Δ, dl, dt, ratio, kw...) |> pairs |> OrderedDict

end
update = update
setup = setup