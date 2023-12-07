using CairoMakie
using CairoMakie: Axis
function saveimg(integrator)
    @unpack t, p, u = integrator
    saveimg(u, p, t)
end
function saveimg(u::AbstractArray, p::AbstractArray, t)

    # function saveimg(u, p, t)
    #     ignore() do
    #         @unpack ϵ, μ, σ, σm = p
    ϵ, μ, σ, σm = eachslice(p, dims=ndims(p))
    Ez, Hx, Hy, Jz = eachslice(u)
    fig = Figure()
    ax = Axis(fig[1, 1]; title="t = $t\nEz")
    heatmap!(ax, Ez, colormap=:seismic, colorrange=(-0.5, 0.5))
    heatmap!(ax, Array(σ), alpha=0.2, colormap=:binary, colorrange=(0, 4),)
    heatmap!(ax, Array(ϵ), alpha=0.5, colormap=:speed, colorrange=(ϵ1, ϵ2))

    a = zeros(size(ϵ))
    for m = monitor_configs
        a[m.idxs...] = 1
    end
    heatmap!(ax, a, alpha=0.6, colormap=:speed, colorrange=(0, 1))
    # fig = heatmap(a.E.z,colorrange=(-.5,.5))
    save("temp/$t.png", fig)
end