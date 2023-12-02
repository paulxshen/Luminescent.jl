using CairoMakie
using CairoMakie: Axis
function saveimg(integrator)
    ignore() do
        @unpack t, u = integrator
        fig = Figure()
        ax = Axis(fig; title="t = $t\nEz")
        heatmap!(ax, u.Ez, colormap=:seismic, colorrange=(-0.5, 0.5))
        heatmap!(ax, Array(σ), alpha=0.2, colormap=:binary, colorrange=(0, 4),)
        heatmap!(ax, Array(ϵ), alpha=0.5, colormap=:speed, colorrange=(ϵ1, ϵ2))

        a = zeros(size(ϵ))
        for m = monitor_configs
            a[m.pos...] .= 1
        end
        heatmap!(ax, a, alpha=0.2, colormap=:speed, colorrange=(0, 1))
        # fig = heatmap(a.E.z,colorrange=(-.5,.5))
        save("temp/$t.png", fig)
    end
end