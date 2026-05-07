function draw_bbox!(ax, bbox)
    isnothing(bbox) && return
    x0, y0 = bbox[:, 1]
    x1, y1 = bbox[:, 2]
    lines!(ax, [(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)], color=:green)
end
function drawgrid!(ax, rulers)
    xmin, xmax = rulers[1][1], rulers[1][end]
    ymin, ymax = rulers[2][1], rulers[2][end]
    kwargs = (color=:gray, linewidth=1 / PIXELS)
    for x = rulers[1]
        lines!(ax, [(x, ymin), (x, ymax)]; kwargs...)
    end
    for y = rulers[2]
        lines!(ax, [(xmin, y), (xmax, y)]; kwargs...)
    end
end

function plotfield(args...; iview, rulers, bbox, labels=[], show_grid=false, show_box=true, t=nothing, slices, title, kwargs...)
    fig = args[1]
    u = args[end]
    N = ndims(u)
    ps = 0

    ps = []
    for (j, (i, v)) in enumerate(slices)
        q = deleteat!(collect(1:3), i)
        rs = rulers[q]
        ax = Axis(fig[iview, j]; title, aspect=(rs[1][end] - rs[1][1]) / (rs[2][end] - rs[2][1]))

        if isnothing(v)
            # n = round(Int, size(u, i) / 2)
        else
            n = round(Int, indexof(rulers[i], v))
            n = clamp(n, 1, size(u, i))
            a = selectdim(u, i, n)
            q = deleteat!(collect(1:3), i)
            rs = rulers[q]
            heatmap!(ax, rs..., a; kwargs...)
            show_box && draw_bbox!(ax, bbox[q, :])
            show_grid && drawgrid!(ax, rs)
            if i == 3
                for (text, pos) in labels
                    text!(ax, pos[1:2]...; text, align=(:center, :center))
                end
            end
        end
    end
end
pow(x, n) = sign(x) * abs(x)^n
function regfy(a, deltas)
    L = sum.(deltas)
    dl = quantile.(deltas, 0.2)
    sz = ceil.(Int, L ./ dl) |> Tuple
    dl = L ./ sz
    rs = cumsum.(deltas)
    # map(CartesianIndices(sz)) do I
    # tmap(collect(CartesianIndices(sz))) do I
    #     I = Tuple(I)
    #     p = (I - 0.5f0) * dl
    #     J = round.(Int, indexof.(rs, p) + 0.5f0)
    #     J = min.(J, size(a))
    #     a[J...]
    # end
    # @show dl, sz, L, rs
    # error("stop")
    axs = map(rs, sz, dl) do r, n, d
        p = ((1:n) - 0.5f0) * d
        indexof.((r,), p) + 0.5f0
    end
    L, resize(a, axs)
end
function nearest(x, v)
    i = searchsortedfirst(v, x)
    i == lastindex(v) + 1 && return v[end]
    i == 1 && return v[1]
    v[i] - x < x - v[i-1] && return v[i]
    v[i-1]
end
function plotmax(u, g; do3d=false, title, rulers, bbox, Ibbox, label_pos, path, views, px_per_unit=4, azimuth=-120°, elevation=60°, boundaries, kwargs...)
    fig = Figure()#resolution=(800, 600))
    u, umax = u
    g, gmin, gmax = g
    kus = []


    # @debug title
    for (iview, v) in enumerate(views)
        @unpack field_color_intensity, material_color_map, material_color_intensity, show_grid, show_box, labels, mirrors = merge(v, kwargs)


        labels = kvmap(labels) do k, v
            string(v) => label_pos(k)
        end |> pairs

        ku = v(:field) |> Symbol
        # k1 = string(v(:field))[1:1]
        slices = sort(v(:slices), by=x -> -x[1])
        a = u(ku) |> Float32#* 100
        N = ndims(a)

        kg = v(:prop) |> togreek
        b = g(kg) |> Float32

        _rulers = deepcopy(rulers)
        for s = mirrors
            dims, side = string(s)
            dims = findfirst(dims, "xyz")
            v = _rulers[dims]
            if side == '-'
                a = cat(reverse(a; dims), a; dims)
                b = cat(reverse(b; dims), b; dims)
                _rulers[dims] = vcat(v[1] - reverse(cumsum(diff(v))), v)
            elseif side == '+'
                a = cat(a, reverse(a; dims); dims)
                b = cat(b, reverse(b; dims); dims)
                _rulers[dims] = vcat(v, cumsum(diff(v)) + v[end])
            else
                error("mirror $s not recognized")
            end
        end

        colormap = Symbol(ku) ∈ (:E, :H) ? [:transparent, :orange] : [:blue, :transparent, :red]
        colorrange = (Symbol(ku) ∈ (:E, :H) ? (0, umax(ku)) : (-umax(ku), umax(ku))) / field_color_intensity
        # @show (; field_color_intensity, colorrange)

        if isnothing(material_color_map)
            material_color_map = cgrad([colorant"transparent", (:black, 0.5)], [gmin(kg), gmax(kg)])
        else
            material_vals = material_color_map[2]
            material_color_map = material_color_map[1], min.(gmax(kg), material_color_map[2])
            material_color_map = material_color_map |> splat(cgrad)
        end
        geometry_color_range = (gmin(kg), gmax(kg)) / material_color_intensity

        subtitle = "$ku $kg"

        if do3d && ku ∉ kus
            a = a[Ibbox...]
            b = b[Ibbox...]

            deltas = getindex.(diff.(_rulers), Ibbox)
            L, a = regfy(a, deltas)
            # ax = Axis3(fig; title=subtitle, aspect=Tuple(L),)# width=800, height=600)
            # GLMakie.volume!(ax, a; colormap, colorrange, algorithm=:absorption, absorption=field_color_intensity)

            ax = Axis3(fig[1, 1]; title=subtitle, aspect=Tuple(L), elevation, azimuth)
            GLMakie.volume!(ax, a; colormap, colorrange, algorithm=:absorption,)
            #  absorption=field_color_intensity)

            L, b = regfy(b, deltas)
            b = min.(b, geometry_color_range[2])
            @debug extrema(b)

            b = nearest.(b, (material_vals,))
            GLMakie.volume!(ax, b; colormap=material_color_map, colorrange=geometry_color_range, algorithm=:absorption, absorption=material_color_intensity, interpolate=false)

            ax.elevation[] = elevation
            ax.azimuth[] = azimuth
            # error("stop")
            push!(kus, ku)

        else
            plotfield(fig, a; rulers=_rulers, bbox, colormap, colorrange, labels, alpha=1, slices, show_grid, show_box, iview, title=subtitle)

            plotfield(fig, b; rulers=_rulers, bbox, colormap=material_color_map, colorrange=geometry_color_range, alpha=0.5, labels, slices, show_grid, show_box, iview, title=subtitle)

            Label(fig[0, :], title)
        end
    end
    resize_to_layout!(fig)
    # resize!(fig, 800, 600)
    if do3d
        GLMakie.save(joinpath(path), fig; px_per_unit)
        GC.gc()
    else
        CairoMakie.save(joinpath(path), fig; px_per_unit)
    end
    fig
end
Base.parse(::Type{Colorant}, t::Union{Tuple,AbstractVector}) = alphacolor(parse(Colorant, t[1]), t[2])
function makemovie(path; full=true, tmax=Inf, kwargs...)
    GRID = joinpath(path, "grid.json")
    PROB = joinpath(path, "problem.json")
    VIS = joinpath(path, "visualization.json")
    if !ispath(VIS)
        println("no views found, skipping movie generation.")
        return
    end

    vis = readjsonnp(VIS)
    prob = readjsonnp(PROB)
    grid = readjsonnp(GRID)
    println("""
    $BREAK
    generating movie...
    """)


    @unpack views, = vis
    @unpack rulers, bbox, Ibbox, L, N, boundaries = grid
    @unpack monitors, sources, name = prob

    for v = views
        v[:prop] = togreek(v(:prop))
    end

    λ = prob(:wavelength) |> FF
    label_pos = map(vcat(monitors, sources)) do x
        (x(:name) => x(:origin))
    end |> Dict
    rulers = rulers * λ
    bbox = bbox * λ
    L *= λ
    if N == 3
        for v = views
            v[:slices] = [string(v) == "mid" ? (i, L[i] / 2 + bbox[i, 1]) : (i, v) for (i, v) = v(:slices)]
        end
    end
    FIELDS = joinpath(path, "fields")
    FRAMES = joinpath(path, "frames")
    FRAMES_SPECIAL = joinpath(path, "frames_special")
    rm(FRAMES, force=true, recursive=true)
    rm(FRAMES_SPECIAL, force=true, recursive=true)
    mkpath(FRAMES)
    mkpath(FRAMES_SPECIAL)

    g = npzread(joinpath(path, "g.npz"))
    gmax = kvmap(g) do k, v
        if Symbol(k) == :ϵ
            b = isPEC(v)
            v = any(b) ? 2maximum(v[.!b]) : maximum(v)
        else
            v = maximum(v)
        end
        k => v
    end
    gmin = vmap(minimum, g)
    # @show gmax



    ss = full ? ("", "_special") : ("_special",)
    for s = ss
        _FIELDS = "$FIELDS$s"
        _FRAMES = "$FRAMES$s"
        run(`rm -rf $_FRAMES`)
        mkpath(_FRAMES)
        v = readdir(_FIELDS, join=true)
        if !isempty(v)

            v = sort(v, by=x -> parse(Float32, basename(x)[1:end-4]))
            ts = [parse(Float32, basename(x)[1:end-4]) for x in v]

            us = npzread.(v)
            for u = us
                for k = ("E", "H")
                    u[k] = Porcupine.norm(u(Regex("$k(.)")))
                end
            end
            us = kmap.(Symbol, us)
            umax = reduce(us; init=vmap(maximum, us[1])) do d1, d2
                NamedTuple([k => max(maximum(d1(k)), maximum(d2(k))) for k = keys(d1)])
            end

            # tmap(collect(enumerate(zip(us, ts)))) do (i, (u, t))
            map(enumerate(zip(us, ts))) do (i, (u, t))
                # if t <= tmax
                title = "$(titlecase(replace(name,"_" => " ")))\nt = $t"
                args = (u, umax), (g, gmin, gmax)
                # local kwargs, v, bbox, rulers
                kwargs = (; kwargs..., rulers, label_pos, Ibbox, bbox, boundaries, path=joinpath(_FRAMES, "$t.png"), title, views)
                plotmax(args...; kwargs...)
                # end
            end
        end
    end

    println("\nmovie frames completed.")
    0

    # v = sort(readdir(FRAMES, join=true), by=x -> parse(Float32, basename(x)[1:end-4]))
    # encoder_options = (crf=23, preset="medium")
    # fps = max(1, round(Int, fps))g
    # open_video_out(MOVIE, load(first(v)) .|> RGB; fps, encoder_options) do writer
    #     for f = v
    #         write(writer, load(f) .|> RGB)
    #     end
    # end
end
