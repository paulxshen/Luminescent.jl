# try
#     using GLMakie
#     using GLMakie: volume
#     global gl = true
# catch e
# using CairoMakie
global gl = false

° = π / 180
function _plot!(g, a, ; colorrange=nothing, title="", colormap=:seismic, algorithm=nothing,
    azimuth=75°, elevation=75°,
    kw...)
    if ndims(a) == 2 || !gl
        if isnothing(colorrange)
            colorrange = 2extrema(a)
        end
        if ndims(a) == 3
            println("3D array: plotting middle slice")
            title *= " (middle slice of 3D array)"
            a = a[:, :, round(Int, size(a, 3) / 2)]
        end
        aspect = size(a, 1) / size(a, 2)
        ax, pl = heatmap(g[1, 1], real(a); axis=(; kw..., title, aspect), colormap, colorrange=colorrange)
    else
        if isnothing(colorrange)
            colorrange = extrema(a) * 0.1
        end
        ax, pl = volume(g[1, 1], real(a), ; axis=(; kw..., type=Axis3, title,), colormap, colorrange, algorithm)
        ax.elevation[] = elevation
        ax.azimuth[] = azimuth
    end
    if diff(collect(extrema(a)))[1] > 0
        Colorbar(g[1, 2], pl)
    end
end
function quickie(sol; kw...)
    @unpack fields, geometry = sol
    fig = Figure()
    fields = (; H=(; Hz=fields.H.Hz))
    geometry = (; ϵ=geometry.ϵ)
    colorrange = (-1, 1) .* maximum(maximum.(a -> abs.(real(a)), flatten(fields)))
    colormap = :seismic
    algorithm = :absorption
    i = 1
    for k1 = keys(fields)
        j = 1
        for (k2, a) = pairs(fields[k1])
            g = fig[i, j]
            title = string(k2)
            _plot!(g, a, ; title, colormap, colorrange, algorithm, kw...)

            j += 1
        end
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