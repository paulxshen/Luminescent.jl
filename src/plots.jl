using CairoMakie
function plotslices(a; saturation=1, path=nothing, n=1)
    fig = Figure()
    v = maximum(abs, a)
    colormap = :seismic
    colorrange = (-v, v) ./ saturation


    for i in 1:3
        for j = 1:n
            u = selectdim(a, i, round(Int, size(a, i) * j / (n + 1)))
            axis = (; aspect=size(u, 1) / size(u, 2))
            heatmap(fig[i, j], u; colormap, colorrange, axis)
        end
    end

    display(fig)
    !isnothing(path) && CairoMakie.save(path, fig)
end

# a = randn(10, 10, 10)
# plotslices(a)