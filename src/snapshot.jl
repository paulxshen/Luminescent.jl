function draw_bbox!(ax, bbox)
    # isnothing(bbox) && return
    # x0, y0 = bbox[:, 1]
    # x1, y1 = bbox[:, 2]
    # lines!(ax, [(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)], color=:black)
end

° = π / 180
function quickie(u, g=nothing; dl=1, monitor_instances=[], source_instances=[], ulims=nothing, λ=1, unit="um", ratio, bbox=zeros(3), origin=zeros(3),)

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
        # text = isempty(m.label) ? "o$i" : m.label
        text = "o$i"
        push!(labels, (ratio * m.center, text))
    end
    for (i, s) = enumerate(source_instances)
        # text = isempty(s.label) ? "s$i" : s.label
        text = "s$i"
        push!(labels, (ratio * s.center, text))
    end

    ax = fig[1, 1]

    # o=origin - bbox[:, 1] * dl
    o = zeros(3)
    xtickformat = x -> string.(Base.round.(x * dl * λ + o[1], digits=2)) .* (unit,)
    ytickformat = y -> string.(Base.round.(y * dl * λ + o[2], digits=2)) .* (unit,)
    ztickformat = z -> string.(Base.round.(z * dl * λ + o[3], digits=2)) .* (unit,)
    if N == 3
        l, w, h = size(u)
    else
        l, w = size(u)

    end

    title0 = "Hz"
    if N == 3
        ax0 = fig[1, 1]
        ax = ax0[1, 1]
        title = "$title0 (xy slice of 3D array)"
        aspect = l / w
        axis = (; title, aspect, xtickformat, ytickformat)
        a = u[:, :, round.(Int, size(u, 3) / 2)]
        heatmap(ax, a; axis, colormap, colorrange)
        draw_bbox!(ax, bbox)

        title = "$title0 (yz slice of 3D array)"
        ax0 = fig[1, 2]
        ax = ax0[1, 1]
        aspect = w / h
        axis = (; title, aspect, xtickformat=ytickformat, ytickformat=ztickformat)
        a = u[round.(Int, size(u, 1) / 2), :, :]
        heatmap(ax, a; axis, colormap, colorrange=colorrange)
        draw_bbox!(ax, bbox[2:3, :])

    else
        ax0 = fig[1, 1]
        ax = ax0[1, 1]
        title = " $title0 (2D array)"
        aspect = l / w
        axis = (; title, aspect, xtickformat, ytickformat)
        heatmap(ax, u; axis, colormap, colorrange)
        draw_bbox!(ax, bbox)

        ax = ax0[1, 2]
        # if !isnothing(g)
        #     if diff(collect(extrema(u)))[1] > 0
        #         Colorbar(ax[1, 2], plt)
        #         contour = g .> 0.99maximum(g)
        #         contour = morpholaplace(contour,)
        #         heatmap!(ax, contour; colormap=[(:gray, 0), :black], colorrange=(0, 1))
        #     end
        # end
    end
    for (pos, text) in labels
        ax = fig[1, 1][1, 1]
        !isnothing(pos) && text!(ax, pos..., ; text, align=(:center, :center))
        # annotate!(g[1, 1], pos, text; fontsize=10, color=:black)
    end
    if isnothing(g)
        return fig
    end

    title0 = "ϵ"
    colormap = [:white, :gray]
    if N == 3
        ax0 = fig[2, 1]
        ax = ax0[1, 1]
        title = "$title0 (xy slice of 3D array)"
        aspect = l / w
        axis = (; title, aspect, xtickformat, ytickformat)
        a = g[:, :, round.(Int, size(g, 3) / 2)]
        heatmap(ax, a; axis, colormap)
        draw_bbox!(ax, bbox)

        title = "$title0 (yz slice of 3D array)"
        ax0 = fig[2, 2]
        ax = ax0[1, 1]
        aspect = w / h
        axis = (; title, aspect, xtickformat=ytickformat, ytickformat=ztickformat)
        a = g[round.(Int, size(g, 1) / 2), :, :]
        heatmap(ax, a; axis, colormap,)
        draw_bbox!(ax, bbox[2:3, :])
    else
        ax0 = fig[2, 1]
        ax = ax0[1, 1]
        title = " $title0 (2D array)"
        aspect = l / w
        axis = (; title, aspect, xtickformat, ytickformat)
        _, plt = heatmap(ax, g; axis, colormap)
        draw_bbox!(ax, bbox)

        ax = ax0[1, 2]
        Colorbar(ax, plt)
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
#             title = " $title0 (middle slice of 3D array)"

#             a = a[:, :, round(Int, size(a, 3) / 2)]
#             ax, plt = heatmap(g[1, 1], real(a); axis=(;  title, aspect), colormap, colorrange=colorrange)

#             a = a[round(Int, size(a, 1) / 2), :, :]
#             ax, plt = heatmap(g[1, 2], real(a); axis=(;  title, aspect), colormap, colorrange=colorrange)
#         else
#             title = " $title0 (2D array)"
#             ax, plt = heatmap(g[1, 1], real(a); axis=(;  title, aspect), colormap, colorrange=colorrange)
#         end
#         for (pos, text) in labels
#             text!(g[1, 1], pos..., ; text, align=(:center, :center))
#             # annotate!(g[1, 1], pos, text; fontsize=10, color=:black)
#         end
#     else
#         if isnothing(colorrange)
#             colorrange = extrema(a) * 0.1
#         end
#         ax, plt = volume(g[1, 1], real(a), ; axis=(;  type=Axis3, title,), colormap, colorrange, algorithm)
#         ax.elevation[] = elevation
#         ax.azimuth[] = azimuth
#     end
#     if diff(collect(extrema(a)))[1] > 0
#         Colorbar(g[1, d], plt)
#     end
# end
