function saveimg(integrator)
    ignore() do
        @unpack t, u = integrator
        fig = heatmap(u.E.z,
            axis=(; title="t = $t\nEz"),
            colormap=:seismic, colorrange=(-0.5, 0.5))
        heatmap!(ϵ, alpha=0.5, colormap=:speed, colorrange=(ϵ1, ϵ2))
        heatmap!(σ, alpha=0.2, colormap=:binary, colorrange=(0, 4))

        a = zeros(size(ϵ))
        for m = monitor_configs
            a[m.pos...] .= 1
        end
        heatmap!(a, alpha=0.2, colormap=:speed, colorrange=(0, 1))
        # fig = heatmap(a.E.z,colorrange=(-.5,.5))
        save("temp/$t.png", fig)
    end
end