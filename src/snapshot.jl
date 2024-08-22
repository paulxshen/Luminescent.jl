# try
#     using GLMakie
#     using GLMakie: volume
#     global gl = true
# catch e
# using CairoMakie
global gl = false

° = π / 180
function _plot!(g, a, ; colorrange=nothing, title="", labels=[], colormap=:seismic, algorithm=nothing,
    azimuth=75°, elevation=75°,
    kw...)
    d = ndims(a)
    if d == 2 || !gl
        if isnothing(colorrange)
            colorrange = extrema(a)
        end
        aspect = size(a, 1) / size(a, 2)
        if d == 3
            println("3D array: plotting middle slice")
            title *= " (middle slice of 3D array)"

            a1 = a[:, :, round(Int, size(a, 3) / 2)]
            ax, pl = heatmap(g[1, 1], real(a1); axis=(; kw..., title, aspect), colormap, colorrange=colorrange)

            a1 = a[round(Int, size(a, 1) / 2), :, :]
            ax, pl = heatmap(g[1, 2], real(a1); axis=(; kw..., title, aspect), colormap, colorrange=colorrange)
        else
            title *= " (2D array)"
            ax, pl = heatmap(g[1, 1], real(a); axis=(; kw..., title, aspect), colormap, colorrange=colorrange)
        end
        for (pos, text) in labels
            text!(g[1, 1], pos..., ; text, align=(:center, :center))
            # annotate!(g[1, 1], pos, text; fontsize=10, color=:black)
        end
    else
        if isnothing(colorrange)
            colorrange = extrema(a) * 0.1
        end
        ax, pl = volume(g[1, 1], real(a), ; axis=(; kw..., type=Axis3, title,), colormap, colorrange, algorithm)
        ax.elevation[] = elevation
        ax.azimuth[] = azimuth
    end
    if diff(collect(extrema(a)))[1] > 0
        Colorbar(g[1, d], pl)
    end
end
function quickie(sol, prob; kw...)
    @unpack fields, geometry = sol |> cpu
    @unpack monitor_instances, source_instances = prob |> cpu
    fig = Figure()
    fields = (; Hz=fields.Hz)
    geometry = (; ϵ=geometry.ϵ)
    colorrange = (-1, 1) .* maximum(maximum.(a -> abs.(real(a)), leaves(fields)))
    colormap = :seismic
    algorithm = :absorption
    labels = []
    for (i, m) = enumerate(monitor_instances)
        text = isempty(m.label) ? "o$i" : m.label
        push!(labels, (m.center, text))
    end
    for (i, s) = enumerate(source_instances)
        text = isempty(s.label) ? "s$i" : s.label
        push!(labels, (s.center, text))
    end

    i = 1
    for k = keys(fields)
        j = 1
        # for (k2, a) = pairs(fields[k1])
        a = fields[k]
        g = fig[i, j]
        title = string(k)
        _plot!(g, a, ; title, colormap, colorrange, algorithm, labels, kw...)

        # j += 1
        i += 1
    end
    if !isnothing(geometry)
        j = 1
        for (k2, a) = pairs(geometry)
            g = fig[i, j]
            title = string(k2)
            _plot!(g, a, ; title, algorithm=:mip, kw...)

            j += 1
        end
    end
    return fig
end