function calibrate_mode(mode, ϵmode2d::AbstractVector, dx)
    F = Float32

    ϵ = mean(ϵmode2d)
    T = 8
    l = T / ϵ
    T1 = T + 2
    Δ = 2
    T2 = T1 + Δ
    # wm = wwg / 2
    # w = 2wm + wwg
    # sz = (round(Int, l / dx), round(Int, w / dx))
    ϵmin = minimum(ϵmode2d)
    ws = length(ϵmode2d) * dx
    w = 1.5ws
    sz = round.(Int, (l / dx, w / dx))

    μ = ones(F, sz)
    σ = zeros(F, sz)
    σm = zeros(F, sz)

    ϵ = repeat(ϵmode2d, 1, round(Int, l / dx))'
    top = round((sz[2] - length(ϵmode2d)) / 2)
    bottom = sz[2] - top - length(ϵmode2d)
    ϵ = pad(ϵ, :replicate, (0, bottom), (0, top))
    # ϵ[1+round(wm / dx):end-round(wm / dx)] .= ϵ2

    # wm = 0.8w
    m = Monitor([l, w / 2], [0, ws], [1, 0],)
    # s = Source([0, w / 2], [0, wm]; EH2JM(mode)...)
    s = ModalSource(t -> cispi(2t), [0.1, w / 2], [0, -ws / 2], [0, ws / 2]; EH2JM(mode)...)
    configs = maxwell_setup([], [s], [m], dx, sz; F, ϵmin)
    @unpack dx, dt, sz, geometry_padding, geometry_staggering, field_padding, source_instances, monitor_instances, u0, field_names = configs
    global configs

    A = support(monitor_instances[1])
    p = apply(geometry_padding; ϵ, μ, σ, σm)
    p = apply(geometry_staggering, p)

    u = reduce((u, t) -> maxwell_update!(u, p, t, dx, dt, field_padding, source_instances;), 0:dt:T1, init=deepcopy(u0))
    # fp0 = zeros(F, length(fields(s)))
    u, fp, pp = reduce(T1+dt:dt:T2, init=(u, 0, 0)) do (u, fp, pp), t
        (maxwell_update(u, p, t, dx, dt, field_padding, source_instances),
            fp + cispi(2t) * sqrt(dt / Δ) * dict([k => field(u, k, monitor_instances[1],) for k = field_names]),
            pp + dt / Δ * flux(u, monitor_instances[1],),
        )
    end

    # fp = dict([k => v for (k, v) = zip(field_names, fp)])
    mode = (Ex=fp.Ey |> reverse, Ez=fp.Ex |> reverse)
    port_powers = mean(pp) .* A
    mode, port_powers, pp
end

function collapse_mode(m, ϵ)
    @unpack Ex, Ey, Ez, Hx, Hy, Hz = m
    E = sqrt.(abs2.(Ex) + abs2.(Ey) + abs2.(Ez))
    Ex, Ez, Hy = map([Ex, Ez, Hy]) do a
        mean(a, dims=2) |> vec
    end
    (; Ex, Ez), sum(E .* ϵ, dims=2) ./ sum(E, dims=2) |> vec
end