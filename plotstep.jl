using CairoMakie
using CairoMakie: Axis
function plotstep(integrator)
    @unpack t, p, u = integrator
    plotstep(u, p, t)
end
function plotstep(u::AbstractArray, p::AbstractArray, t; ax=nothing, colorrange)

    # ϵ, μ, σ, σm = p
    # Ez, Hx, Hy, Jz = u
    ϵ, μ, σ, σm = eachslice(p, dims=ndims(p))
    Ez, Hx, Hy, Jz = eachslice(u, dims=ndims(u))
    a = deepcopy(Array(Ez))
    heatmap!(ax, a; colormap=:seismic, colorrange)
    # heatmap!(ax, Array(σ), alpha=0.2, colormap=:binary, colorrange=(0, 4),)
    heatmap!(ax, Array(ϵ), colormap=[(:white, 0), (:gray, 0.6)], colorrange=(ϵ1, ϵ2))

    marker = :rect
    for m = monitor_configs
        x, y = m.idxs
        scatter!(ax, [x], [y]; marker)
    end
    scatter!(ax, [1], [1])
    # heatmap!(ax, a, alpha=0.6, colormap=:speed, colorrange=(0, 1))
    # fig = heatmap(a.E.z,colorrange=(-.5,.5))
    # save("temp/$t.png", fig)
end