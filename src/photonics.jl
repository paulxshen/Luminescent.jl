function calibrate_mode(mode, ϵmode, dx, λ=1.55; F=Float32, verbose=false, name="", kwargs...)
    mode, ϵmode, dx, λ = (mode, ϵmode, dx, λ) |> F
    d = ndims(ϵmode) + 1
    n = sqrt(mean(ϵmode))

    source_monitor_margin, source_boundary_margin, mode_margin = whole.((SOURCE_MONITOR_MARGIN, SOURCE_BOUNDARY_MARGIN, MODE_MARGIN), dx) / λ
    l0 = 0.8
    l = l0 + source_boundary_margin + source_monitor_margin

    # wm = wwg / 2
    # w = 2wm + wwg
    # sz = (round(Int, l / dx), round(Int, w / dx))
    L = [size(ϵmode) * dx...]
    sz = (round(l / dx), size(ϵmode)...)

    ϵ = stack(fill(ϵmode, sz[1]))
    if d == 2
        ϵ = ϵ'
    else
        ϵ = permutedims(ϵ, (3, 1, 2))
    end
    # n = round(mode_margin / dx)
    # ϵ = pad(ϵ, :replicate, (0, n), (0, n))
    # ϵ |> heatmap |> display
    # ϵ[1+round(wm / dx):end-round(wm / dx)] .= ϵ2

    normal = [1, zeros(d - 1)...]
    tangent = [0, -1, zeros(d - 2)...]
    m = [ModalMonitor(mode, [x, (L / 2)...], normal, tangent, L) for x = range(0, l0, 3) + source_monitor_margin + source_boundary_margin]
    s = ModalSource(t -> cispi(2t), mode, [source_monitor_margin, (L / 2)...], -normal, tangent, L;)
    prob = setup([], [s], m, dx, sz; F, ϵ, verbose)
    global aa = prob
    sol = solve(prob; comprehensive=true, verbose=true,)
    @unpack fields, modes, forward_mode_powers = sol
    mode = modes[1][1]
    # global a, b = modes, mode
    power = forward_mode_powers[1][1][1]
    (; mode, power, sol)
    # recordsim("$(@__DIR__)/$name.mp4", h, ;
    #     dt,
    #     field=:Hz,
    #     monitor_instances,
    #     source_instances,
    #     geometry=p.ϵ.ϵyy,
    # )

end
