
° = π / 180
function quickie(u, g=nothing; dx=1, monitor_instances=[], source_instances=[], ulims=nothing, λ=1, unit="um", kw...)
    ratio = dx / 0.025
    u = imresize(u; ratio)

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
        push!(labels, (ratio * m.center, text))
    end
    for (i, s) = enumerate(source_instances)
        text = isempty(s.label) ? "s$i" : s.label
        push!(labels, (ratio * s.center, text))
    end

    grid = fig[1, 1]
    f = x -> string.(round.(x * dx * λ / ratio, digits=2)) .* (unit,)
    ytickformat = xtickformat = f
    if N == 3
        l, w, h = size(u)
    else
        l, w = size(u)

    end
    axis0 = (; kw..., xtickformat, ytickformat)

    title = "Hz"
    if N == 3
        # println("3D array: plotting middle slice")
        title *= " (xy slice of 3D array)"
        aspect = l / w
        axis = merge(axis0, (; title, aspect))
        u = u[:, :, round.(Int, size(u, 3) / 2)]
        ax, plt = heatmap(grid[1, 1], real(u); axis, colormap, colorrange)

        aspect = w / h
        axis = merge(axis0, (; title, aspect))
        u = u[round.(Int, size(u, 1) / 2), :, :]
        ax, plt = heatmap(grid[1, 2], real(u); axis, colormap, colorrange=colorrange)
    else
        title *= " (2D array)"
        aspect = l / w
        axis = merge(axis0, (; title, aspect))
        ax, plt = heatmap(grid[1, 1], u; axis, colormap, colorrange)
        if !isnothing(g)
            if diff(collect(extrema(u)))[1] > 0
                Colorbar(grid[1, 2], plt)
                contour = g .> 0.99maximum(g)
                contour = morpholaplace(contour,)
                heatmap!(ax, contour; colormap=[(:gray, 0), :black], colorrange=(0, 1))
            end
        end
    end
    for (pos, text) in labels
        text!(grid[1, 1], pos..., ; text, align=(:center, :center))
        # annotate!(g[1, 1], pos, text; fontsize=10, color=:black)
    end
    if isnothing(g)
        return fig
    end

    g = imresize(g, size(u))
    grid = fig[2, 1]
    title = "ϵ"
    colormap = [:white, :gray]
    if N == 3
        # println("3D array: plotting middle slice")
        title *= " (middle slice of 3D array)"

        aspect = l / w
        axis = merge(axis0, (; title, aspect))
        g = g[:, :, round.(Int, size(a, 3) / 2)]
        ax, plt = heatmap(grid[1, 1], g; axis, colormap)

        aspect = w / h
        axis = merge(axis0, (; title, aspect))
        g = g[round.(Int, size(a, 1) / 2), :, :]
        ax, plt = heatmap(grid[1, 2], g; axis, colormap,)
    else
        title *= " (2D array)"
        aspect = l / w
        axis = merge(axis0, (; title, aspect))
        ax, plt = heatmap(grid[1, 1], g; axis, colormap)
        if diff(collect(extrema(g)))[1] > 0
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

# function _plot!(g, a, ; colorrange=nothing, title="", labels=[], colormap=:seismic, algorithm=nothing,
#     azimuth=75°, elevation=75°,
#     kw...)
#     d = ndims(a)
#     if d == 2 || !gl
#         if isnothing(colorrange)
#             colorrange = extrema(a)
#         end
#         aspect = size(a, 1) / size(a, 2)
#         if d == 3
#             # println("3D array: plotting middle slice")
#             title *= " (middle slice of 3D array)"

#             a = a[:, :, round(Int, size(a, 3) / 2)]
#             ax, plt = heatmap(g[1, 1], real(a); axis=(; kw..., title, aspect), colormap, colorrange=colorrange)

#             a = a[round(Int, size(a, 1) / 2), :, :]
#             ax, plt = heatmap(g[1, 2], real(a); axis=(; kw..., title, aspect), colormap, colorrange=colorrange)
#         else
#             title *= " (2D array)"
#             ax, plt = heatmap(g[1, 1], real(a); axis=(; kw..., title, aspect), colormap, colorrange=colorrange)
#         end
#         for (pos, text) in labels
#             text!(g[1, 1], pos..., ; text, align=(:center, :center))
#             # annotate!(g[1, 1], pos, text; fontsize=10, color=:black)
#         end
#     else
#         if isnothing(colorrange)
#             colorrange = extrema(a) * 0.1
#         end
#         ax, plt = volume(g[1, 1], real(a), ; axis=(; kw..., type=Axis3, title,), colormap, colorrange, algorithm)
#         ax.elevation[] = elevation
#         ax.azimuth[] = azimuth
#     end
#     if diff(collect(extrema(a)))[1] > 0
#         Colorbar(g[1, d], plt)
#     end
# end
