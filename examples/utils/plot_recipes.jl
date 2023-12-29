using CairoMakie
using CairoMakie: Axis

function plotstep!(ax, u::AbstractVector, p::AbstractVector, configs; colorrange)
    "plot field"
    heatmap!(ax, u[1]; colormap=:seismic, colorrange)
    "plot geometry"
    heatmap!(ax, p[1], colormap=[(:white, 0), (:gray, 0.8)])#, colorrange=(ϵ1, ϵ2))

    "plot monitors"
    marker = :rect
    for (i, m) = enumerate(configs.monitor_info)
        x, y = m.idxs
        scatter!.(ax, x, y, ; marker, color=:yellow)
        x, y = m.center
        text!(ax, x, y, ; text="$i", align=(:center, :center))
    end
    # save("temp/$t.png", fig)
end

function recordsim(sol, p, fdtd_configs, fn; title="", frameat=1 / 16, framerate=16)
    @unpack dt, T = fdtd_configs
    t = 0:frameat:T

    fig = Figure()
    umax = maximum(map(sol) do v
        maximum(abs, v[1])
    end)
    r = record(fig, fn, t; framerate) do t
        i = round(Int, t / dt + 1)
        empty!(fig)
        ax = Axis(fig[1, 1]; title="t = $t\nEz\n$title")
        u = sol[i]
        colorrange = (-1, 1) .* umax
        plotstep!(ax, u, p, fdtd_configs; colorrange)
    end
    println("saved simulation recording to $fn")
    r
end

# function plotmonitors(sol, monitor_info)

#     t = range(0, T, length=size(sol)[end])
#     fig = Figure()
#     for i = 1:2
#         ax = Axis(fig[i, 1], title="Monitor $i")
#         E, H = get(sol, monitor_info[i],)
#         lines!(ax, t, E)
#         lines!(ax, t, H)
#         # lines!(ax, t, H .* E)
#     end
#     display(fig)
# end