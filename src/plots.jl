function plotslices(a; saturation=10, path=nothing)
    fig = Figure()
    v = maximum(abs, a)
    colormap = :seismic
    colorrange = (-v, v) / saturation

    u = a[round(Int, end / 2), :, :]
    axis = (; aspect=size(u, 2) / size(u, 1))
    heatmap(fig[1, 1], u; colormap, colorrange, axis)

    u = a[:, round(Int, end / 2), :]
    axis = (; aspect=size(u, 2) / size(u, 1))
    heatmap(fig[2, 1], u; axis, colormap, colorrange)

    u = a[:, :, round(Int, end / 2)]
    axis = (; aspect=size(u, 2) / size(u, 1))
    heatmap(fig[3, 1], u; colormap, colorrange, axis)

    display(fig)
    !isnothing(path) && CairoMakie.save(path, fig)
end

a = randn(10, 10, 10)
plotslices(a)