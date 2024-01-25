function plotstep!(f, u::AbstractArray, p, configs, monitor_instances=[]; s=1, umax=maximum(abs, u), bipolar=!all(u .>= 0), title="")
    colormap = bipolar ? :seismic : [(:white, 0), (:orange, 1)]
    colorrange = bipolar ? (-1, 1) .* umax : (0, umax)
    algorithm = bipolar ? :absorption : :mip

    d = ndims(u)
    "plot field"
    if d == 1
        #   heatmap!
    elseif d == 2
        pl! = heatmap!
        pl = heatmap(f, u; axis=(; title, aspect=1), colormap, colorrange)
    else
        pl! = volume!
        volume(f, u; axis=(; type=Axis3, title), algorithm, colormap, colorrange)
    end
    # Colorbar(f)

    if !isnothing(p)
        "plot geometry"
        pl!(f, p, colormap=[(:white, 0), (:gray, 0.2)])#, colorrange=(ϵ1, ϵ2))
    end
    "plot monitors"
    if !isempty(monitor_instances)
        a = zeros(size(u))
        for (i, m) = enumerate(configs.monitor_instances)
            a[first(m.idxs)...] .= 1
            text!(f, first(m.centers)..., ; text="$i", align=(:center, :center))
        end
        pl!(f, a, colormap=[(:white, 0), (:yellow, 0.2)])#, colorrange=(ϵ1, ϵ2))
        # # save("temp/$t.png", fig)
    end
end

function recordsim(sol, p, fdtd_configs, fn, monitor_instances=[];
    s=1,
    umax=s * maximum(abs, sol[round(Int, length(sol) / 2)]),
    bipolar=true,
    title="",
    playback=1, frameat=fdtd_configs.dt,
    framerate=playback / frameat)
    @unpack dt, T = fdtd_configs
    t = 0:frameat:T

    fig = Figure()

    r = record(fig, fn, t; framerate) do t
        i = round(Int, t / dt + 1)
        empty!(fig)
        # ax = Axis(fig[1, 1];)
        u = sol[i]

        plotstep!(fig[1, 1], u, p, fdtd_configs, monitor_instances; bipolar, umax, title="t = $t\n$title")
    end
    println("saved simulation recording to $fn")
    r
end
