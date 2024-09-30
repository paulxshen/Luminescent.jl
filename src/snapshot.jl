
° = π / 180
function quickie(u, g=nothing; monitor_instances=[], source_instances=[], ulims=nothing, kw...)
    fig = Figure()
    N = ndims(u)
    if ulims == nothing
        colorrange = (-1, 1) .* maximum(abs, u)
    else
        colorrange = ulims
    end
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

    a = u
    grid = fig[1, 1]
    aspect = size(a, 1) / size(a, 2)
    title = "Hz"
    if N == 3
        # println("3D array: plotting middle slice")
        title *= " (middle slice of 3D array)"

        a1 = a[:, :, round(Int, size(a, 3) / 2)]
        ax, plt = heatmap(grid[1, 1], real(a1); axis=(; kw..., title, aspect), colormap, colorrange=colorrange)

        a1 = a[round(Int, size(a, 1) / 2), :, :]
        ax, plt = heatmap(grid[1, 2], real(a1); axis=(; kw..., title, aspect), colormap, colorrange=colorrange)
    else
        title *= " (2D array)"
        axis = (; kw..., title, aspect)
        ax, plt = heatmap(grid[1, 1], a; axis, colormap, colorrange)
        if diff(collect(extrema(a)))[1] > 0
            Colorbar(grid[1, 2], plt)
            contour = g .> 0.99maximum(g)
            contour = morpholaplace(contour,)
            heatmap!(ax, contour; colormap=[(:gray, 0), :black], colorrange=(0, 1))
        end
    end
    for (pos, text) in labels
        text!(grid[1, 1], pos..., ; text, align=(:center, :center))
        # annotate!(g[1, 1], pos, text; fontsize=10, color=:black)
    end
    if isnothing(g)
        return fig
    end

    a = g
    grid = fig[2, 1]
    aspect = size(a, 1) / size(a, 2)
    title = "ϵ"
    colormap = [:white, :gray]
    if N == 3
        # println("3D array: plotting middle slice")
        title *= " (middle slice of 3D array)"

        a1 = a[:, :, round(Int, size(a, 3) / 2)]
        ax, plt = heatmap(grid[1, 1], real(a1); axis=(; kw..., title, aspect), colormap)

        a1 = a[round(Int, size(a, 1) / 2), :, :]
        ax, plt = heatmap(grid[1, 2], real(a1); axis=(; kw..., title, aspect), colormap,)
    else
        title *= " (2D array)"
        axis = (; kw..., title, aspect)
        ax, plt = heatmap(grid[1, 1], a; axis, colormap)
        if diff(collect(extrema(a)))[1] > 0
            Colorbar(grid[1, 2], plt)
        end
    end

    # i = 1
    # for k = keys(fields)
    #     j = 1
    #     # for (k2, a) = pairs(fields[k1])
    #     a = fields[k]
    #     g = fig[i, j]
    #     title = string(k)
    #     _plot!(g, a, ; title, colormap, colorrange, algorithm, labels, kw...)

    #     # j += 1
    #     i += 1
    # end
    # if !isnothing(geometry)
    #     j = 1
    #     for (k2, a) = pairs(geometry)
    #         g = fig[i, j]
    #         title = string(k2)
    #         _plot!(g, a, ; title, algorithm=:mip, kw...)

    #         j += 1
    #     end
    # end
    return fig
end

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
            # println("3D array: plotting middle slice")
            title *= " (middle slice of 3D array)"

            a1 = a[:, :, round(Int, size(a, 3) / 2)]
            ax, plt = heatmap(g[1, 1], real(a1); axis=(; kw..., title, aspect), colormap, colorrange=colorrange)

            a1 = a[round(Int, size(a, 1) / 2), :, :]
            ax, plt = heatmap(g[1, 2], real(a1); axis=(; kw..., title, aspect), colormap, colorrange=colorrange)
        else
            title *= " (2D array)"
            ax, plt = heatmap(g[1, 1], real(a); axis=(; kw..., title, aspect), colormap, colorrange=colorrange)
        end
        for (pos, text) in labels
            text!(g[1, 1], pos..., ; text, align=(:center, :center))
            # annotate!(g[1, 1], pos, text; fontsize=10, color=:black)
        end
    else
        if isnothing(colorrange)
            colorrange = extrema(a) * 0.1
        end
        ax, plt = volume(g[1, 1], real(a), ; axis=(; kw..., type=Axis3, title,), colormap, colorrange, algorithm)
        ax.elevation[] = elevation
        ax.azimuth[] = azimuth
    end
    if diff(collect(extrema(a)))[1] > 0
        Colorbar(g[1, d], plt)
    end
end
