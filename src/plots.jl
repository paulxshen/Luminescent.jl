function plotslices(a, path=nothing)
    fig = Figure()
    v = maximum(abs, a)
    colormap = :seismic
    colorrange = (-v, v)
    heatmap(fig[1, 1], a[round(Int, end / 2), :, :]; colormap, colorrange)
    heatmap(fig[2, 1], a[:, round(Int, end / 2), :]; colormap, colorrange)
    heatmap(fig[3, 1], a[:, :, round(Int, end / 2)]; colormap, colorrange)
    display(fig)
    !isnothing(path) && CairoMakie.save(path, fig)
end

a = randn(10, 10, 10)
plotslices(a, "a.png")