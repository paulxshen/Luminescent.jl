

function plotstep!(f, u::AbstractArray, p, configs; s=1, umax=maximum(abs, u), bipolar=!all(u .>= 0), title="")
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
    if !isempty(configs.monitor_instances)
        a = zeros(size(u))
        for (i, m) = enumerate(configs.monitor_instances)
            a[m.idxs...] .= 1
            text!(f, m.center..., ; text="$i", align=(:center, :center))
        end
        pl!(f, a, colormap=[(:white, 0), (:yellow, 0.2)])#, colorrange=(ϵ1, ϵ2))
        # # save("temp/$t.png", fig)
    end
end

function recordsim(sol, p, fdtd_configs, fn;
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

        plotstep!(fig[1, 1], u, p, fdtd_configs; bipolar, umax, title="t = $t\n$title")
    end
    println("saved simulation recording to $fn")
    r
end

# function plotmonitors(sol, monitor_instances)

#     t = range(0, T, length=size(sol)[end])
#     fig = Figure()
#     for i = 1:2
#         ax = Axis(fig[i, 1], title="Monitor $i")
#         E, H = get(sol, monitor_instances[i],)
#         lines!(ax, t, E)
#         lines!(ax, t, H)
#         # lines!(ax, t, H .* E)
#     end
#     display(fig)
# end