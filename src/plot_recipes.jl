using CairoMakie
using CairoMakie: Axis
function plotstep(integrator)
    @unpack t, p, u = integrator
    plotstep(u, p, t)
end
function plotstep(u::AbstractArray, p::AbstractArray, t, configs; ax=nothing, colorrange)

    ϵ, μ, σ, σm = p
    Ez, Hx, Hy = u
    # ϵ, μ, σ, σm = eachslice(p, dims=ndims(p))
    # Ez, Hx, Hy, Jz = eachslice(u, dims=ndims(u))
    a = deepcopy(Array(Ez))
    heatmap!(ax, a; colormap=:seismic, colorrange)
    # heatmap!(ax, Array(σ), alpha=0.2, colormap=:binary, colorrange=(0, 4),)
    heatmap!(ax, Array(ϵ), colormap=[(:white, 0), (:gray, 0.6)])#, colorrange=(ϵ1, ϵ2))

    marker = :rect
    for m = configs.monitor_configs
        x, y = m.idxs
        scatter!(ax, [x], [y]; marker)
    end
    scatter!(ax, [1], [1])
    # heatmap!(ax, a, alpha=0.6, colormap=:speed, colorrange=(0, 1))
    # fig = heatmap(a.E.z,colorrange=(-.5,.5))
    # save("temp/$t.png", fig)
end

function recordsim(sol, p, fdtd_configs, fn; title="", frameat=1 / 16, framerate=16)
    @unpack dt, T = fdtd_configs


    t = 0:frameat:T

    fig = Figure()
    umax = maximum(sol[end][1])
    record(fig, fn, t; framerate) do t
        # i = round.(Int, t ./ dt) .+ 1
        i = round(Int, t / dt + 1)
        empty!(fig)
        ax = Axis(fig[1, 1]; title="t = $t\nEz\n$title")
        u = sol[i]
        colorrange = (-1, 1) .* umax
        # u = sol[:, :, :, i]
        plotstep(u, p, t, fdtd_configs; ax, colorrange)
    end
end

using CairoMakie
using CairoMakie: Axis
function plotmonitors(sol, monitor_configs)

    t = range(0, T, length=size(sol)[end])
    fig = Figure()
    for i = 1:2
        ax = Axis(fig[i, 1], title="Monitor $i")
        E, H = get(sol, monitor_configs[i],)
        lines!(ax, t, E)
        lines!(ax, t, H)
        # lines!(ax, t, H .* E)
    end
    display(fig)
end