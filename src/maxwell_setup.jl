function add_current_keys(d)
    # d=OrderedDict
    d = dict(pairs(d))
    add_current_keys!(d)
end
function add_current_keys!(d::AbstractDict)
    for k in keys(d) |> collect
        if startswith(String(k), "E")
            d[Symbol("J" * String(k)[2:end])] = d[k]
        end

        if startswith(String(k), "H")
            d[Symbol("M" * String(k)[2:end])] = d[k]
        end

    end
    d
end
"""
    function maxwell_setup(boundaries, sources, monitors, L, dx, polarization=nothing; F=Float32)

Args
...
- L: vector of lengths in wavelengths of simulation domain
- polarization: only applies to 2d which can be :TM (Ez, Hx, Hy) or :TE (Hz, Ex, Ey)
"""
function maxwell_setup(boundaries, sources, monitors, dx, sz;
    polarization=:TE,
    geometry=OrderedDict(),
    # transient_duration=max_source_dist(sources), steady_state_duration=2,
    transient_duration=nothing, steady_state_duration=2,
    ϵ=1, μ=1, σ=0, σm=0, F=Float32,
    wavelengths=[1],
    # ϵmin=1,
    # Courant=0.5,
    Courant=nothing,
    verbose=true,
    kw...)
    dx, wavelengths = F.((dx, wavelengths))

    a = ones(F, sz)
    geometry = dict((; ϵ, μ, σ, σm) |> pairs)
    for (k, v) = pairs(geometry)
        if isa(v, Number)
            geometry[k] = a * geometry[k]
        end
    end
    ϵmin, ϵmax = extrema(geometry.ϵ)
    if isnothing(Courant)
        Courant = F(0.8√(ϵmin / length(sz)))# Courant number)
    end
    if isnothing(transient_duration)
        transient_duration = sum(sz) * dx * sqrt(ϵmax) + 1
    end

    d = length(sz)
    if d == 1
        field_names = (:Ez, :Hy)
        polarization = nothing
    elseif d == 2
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
    Courant = F(Courant)

    nodes = fill(:U, d, 2)
    db = Any[PML(j * i,) for i = 1:d, j = (-1, 1)]
    field_padding = DefaultDict(() -> InPad[])
    geometry_padding = DefaultDict(() -> OutPad[])
    fl = Dict([k => zeros(Int, d) for k = field_names])
    is_field_on_lb = Dict([k => zeros(Int, d) for k = field_names])
    is_field_on_ub = Dict([k => zeros(Int, d) for k = field_names])
    common_left_pad_amount = zeros(Int, d)
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

    for i = 1:d
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

    for i = 1:d
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



    for i = 1:d
        select = i .== 1:d
        xyz = para = perp = [:x, :y, :z]

        perp = [popat!(para, i)]
        for j = 1:2
            b = db[i, j]
            if isa(b, PML)
                pml_amount = round(Int, b.d / dx)
                l = j == 1 ? [i == a ? pml_amount : 0 for a = 1:d] : zeros(Int, d)
                r = j == 2 ? [i == a ? pml_amount : 0 for a = 1:d] : zeros(Int, d)
                push!(geometry_padding[:ϵ], OutPad(:replicate, l, r, sz))
                push!(geometry_padding[:μ], OutPad(:replicate, l, r, sz))

                l1 = round(l / 2)
                r1 = round(r / 2)
                push!(geometry_padding[:σ], OutPad(ReplicateRamp(F(b.σ)), l1, r1, sz))
                push!(geometry_padding[:σm], OutPad(ReplicateRamp(F(b.σ)), l1, r1, sz))
                l2 = l - l1
                r2 = r - r1
                push!(geometry_padding[:σ], OutPad(F(b.σ), l2, r2, sz))
                push!(geometry_padding[:σm], OutPad(F(b.σ), l2, r2, sz))
                if j == 1
                    for k = keys(fl)
                        fl[k][i] += pml_amount
                    end
                    common_left_pad_amount[i] += pml_amount
                end
                for k = keys(field_sizes)
                    field_sizes[k][i] += pml_amount
                end
            end
            l = j == 1 ? Int.((1:d) .== i) : zeros(Int, d)
            r = j == 2 ? Int.((1:d) .== i) : zeros(Int, d)
            t = typeof(b)


            f = nodes[i, j]
            for k = field_names
                q = startswith(String(k), String(f))
                if (q ? k[2] in para : k[2] in perp)
                    if t == Periodic
                        lp = j == 1 ? select : zeros(Int, d)
                        rp = j == 2 ? select : zeros(Int, d)
                        push!(field_padding[k], InPad(:periodic, lp, rp,))
                    else
                        push!(field_padding[k], InPad(0, l, r,))
                    end
                end
            end
        end
    end


    for (k, v) = pairs(field_padding)
        for p = v
            fl[k] += p.l
            is_field_on_lb[k] += p.l
            is_field_on_ub[k] += p.r
        end
    end

    add_current_keys!(fl)
    add_current_keys!(field_sizes)

    field_sizes = NamedTuple([k => Tuple(field_sizes[k]) for (k) = keys(field_sizes)])
    fields = NamedTuple([k => zeros(F, Tuple(field_sizes[k])) for (k) = field_names])
    if d == 1
        pf = u -> [u[1] .* u[2]]
    elseif d == 3
        u0 = dict([
            :E => dict([k => zeros(F, Tuple(field_sizes[k])) for k = (:Ex, :Ey, :Ez)]),
            :H => dict([k => zeros(F, Tuple(field_sizes[k])) for k = (:Hx, :Hy, :Hz)]),
            :J => dict([k => zeros(F, Tuple(field_sizes[k])) for k = (:Jx, :Jy, :Jz)]),
        ])
        # u0 = [u[1:3], u[4:6]]
    else
        if polarization == :TM
            u -> [-u[1] .* u[3], u[2] .* u[1]]
        else
            u0 = dict([
                :E => dict([k => zeros(F, Tuple(field_sizes[k])) for k = (:Ex, :Ey)]),
                :H => dict([k => zeros(F, Tuple(field_sizes[k])) for k = (:Hz,)]),
                :J => dict([k => zeros(F, Tuple(field_sizes[k])) for k = (:Jx, :Jy)]),
            ])
        end
    end

    geometry_sizes = NamedTuple([k => sz .+ sum(geometry_padding[k]) do p
        p.l + p.r
    end for k = keys(geometry_padding)])
    # subpixel_averaging = dict(Pair.(keys(geometry_sizes), Base.oneto.(values(geometry_sizes))))
    subpixel_averaging = Dict{Symbol,Any}()
    field_grids = Dict{Symbol,Any}()
    for k = keys(geometry_sizes)
        if k in (:μ, :σm)
            names = Hnames
        elseif k in (:ϵ, :σ)
            names = Enames
        end
        v = dict([
            begin
                xyz = f[2]
                terminations = zip(is_field_on_lb[f], is_field_on_ub[f])
                g = Symbol("$(k)$xyz$xyz")
                v = SubpixelAveraging(
                    stack([
                        begin
                            if v == (0, 0)
                                [0, 0]
                            elseif v == (0, 1)
                                [0, 0]
                            elseif v == (1, 1)
                                [-1, 1]
                            elseif v == (1, 0)
                                [-1, -1]
                            end
                        end for v in terminations
                    ]))
                if !haskey(field_grids, f)
                    field_grids[f] = v
                end
                g => v
            end for f = names
        ])
        subpixel_averaging[k] = v
    end

    field_origin = NamedTuple([k => (1 - v) * dx / 2 - common_left_pad_amount * dx for (k, v) = pairs(is_field_on_lb)])
    field_origin = add_current_keys(field_origin)

    source_instances = SourceInstance.(sources, dx, (field_sizes,), (field_origin,), (common_left_pad_amount,), (sz,); F)
    monitor_instances = MonitorInstance.(monitors, dx, (field_origin,), (common_left_pad_amount,), (sz,), ; F)
    # roi = MonitorInstance(Monitor(zeros(d), zeros(d), dx * sz), dx, sz, common_left_pad_amount, is_field_on_lb, fl; F)

    dt = dx * Courant
    sz = Tuple(sz)
    for (k, v) = pairs(field_padding)
        a = ones(F, size(fields[k]))
        i = filter(eachindex(v)) do i
            v[i].b == 0
        end
        i_ = sort(setdiff(eachindex(v), i))
        if !isempty(i)

            for i = i
                @unpack b, l, r = v[i]
                pad!(a, b, l, r)
            end
            field_padding[k] = [
                InPad(
                    0,
                    sum(getproperty.(v[i], :l)),
                    sum(getproperty.(v[i], :r)),
                    a),
                v[i_]...
            ]
        end
    end
    geometry_padding = NamedTuple(geometry_padding)
    field_padding = NamedTuple(field_padding)
    if verbose
        @info """
     ====
     FDTD prob
     
     Lengths in characterstic wavelengths, times in characterstic periods, unless otherwise specified

     dx: $(dx|>d2) 
     dt: $(dt|>d2) 
     Courant number: $(Courant|>d2)

     Original array size of all fields in pixels: $sz
     Padded field array field_sizes in pixels:
     $field_sizes
     
     Boundaries:
     $(join("- ".*string.(db),"\n"))
     
     Sources:
     $(join("- ".*string.(sources),"\n"))
     
     Monitors:
     $(join("- ".*string.(monitors),"\n"))
     
     $footer
     ====
     """
    end

    OrderedDict((; geometry_padding, field_padding, subpixel_averaging, field_grids,
        source_instances, monitor_instances, field_names,
        polarization,
        transient_duration, steady_state_duration,
        geometry,
        wavelengths,
        # roi,
        u0, fields, d, dx, dt, sz, kw...) |> pairs)
end

update = maxwell_update
setup = maxwell_setup