function calibrate_mode(mode, ϵmode::AbstractVector, dx; F=Float32, verbose=false, name="")
    mode, ϵmode, dx = (mode, ϵmode, dx) |> F

    n = sqrt(mean(ϵmode))
    l = 1.6
    T1 = l * n + 2
    Δ = 2
    T2 = T1 + Δ
    # wm = wwg / 2
    # w = 2wm + wwg
    # sz = (round(Int, l / dx), round(Int, w / dx))
    ϵmin = minimum(ϵmode)
    ws = length(ϵmode) * dx
    my = round(0.5 / dx) * dx
    w = ws + 2my
    sz = round.(Int, (l / dx, w / dx))

    μ = ones(F, sz)
    σ = zeros(F, sz)
    σm = zeros(F, sz)

    ϵ = repeat(ϵmode, 1, sz[1])'
    top = round(Int, my / dx)
    bottom = sz[2] - top - length(ϵmode)
    ϵ = pad(ϵ, :replicate, (0, bottom), (0, top))
    ϵ |> heatmap |> display
    # ϵ[1+round(wm / dx):end-round(wm / dx)] .= ϵ2

    m = [Monitor([x, w / 2], [0, ws], [1, 0],) for x = 0.5:0.5:1.5]
    s = ModalSource(t -> cispi(2t), [margin_offset, w / 2], [0, ws], [1, 0]; E2J(mode)...)
    configs = maxwell_setup([], [s], m, dx, sz; F, ϵmin, verbose)
    @unpack dx, dt, sz, geometry_padding, geometry_staggering, field_padding, source_instances, monitor_instances, u0, field_names = configs
    # global configs

    A = area(monitor_instances[1])
    p = apply(geometry_padding; ϵ, μ, σ, σm)
    p = apply(geometry_staggering, p)

    u = reduce((u, t) -> maxwell_update!(u, p, t, dx, dt, field_padding, source_instances;), 0:dt:T1, init=deepcopy(u0))
    # fp0 = zeros(F, length(fields(s)))
    u, fp, pp = reduce(T1+dt:dt:T2, init=(u, 0, 0)) do (u, fp, pp), t
        _fp = [invreframe(frame(m), field.((u,), (:Ex, :Ey, :Hz), (m,))) for m = monitor_instances]
        (
            maxwell_update(u, p, t, dx, dt, field_padding, source_instances),
            fp + cispi(-2t) * dt / Δ * _fp,
            pp + dt / Δ * flux.((u,), monitor_instances,),
        )
    end
    fp = [
        begin

            Ex, Hy, Ez = v
            (; Ex, Hy, Ez)
        end for v = fp
    ]

    a = mode_decomp.((mode,), fp)
    a = [
        begin
            ap, am = a
            @show abs(ap)^2, abs(am)^2
        end
        for a = a
    ]
    # recordsim("$(@__DIR__)/$name.mp4", fp.Hz, ;
    #     dt,
    #     field=:Hz,
    #     monitor_instances,
    #     source_instances,
    #     geometry=p.ϵ.ϵEy,
    # )

    @show port_powers = mean.(pp) .* A
    fp[end], a[end][1], pp
end
