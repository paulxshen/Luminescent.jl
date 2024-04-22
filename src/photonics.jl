function calibrate_source(s, dx, wwg, ϵ1, ϵ2)
    F = Float32

    l = 8 / sqrt(ϵ1)
    wm = wwg / 2
    w = 2wm + wwg
    sz = (round(Int, l / dx), round(Int, w / dx))
    ϵmin = ϵ2

    μ = ones(F, sz)
    σ = zeros(F, sz)
    σm = zeros(F, sz)
    ϵ = ones(F, sz) * ϵ1
    ϵ[1+round(wm / dx):end-round(wm / dx)] .= ϵ2

    m = Monitor([l - 0.1, w / 2], [0, -0.75wwg], [0, 0.75wwg]; normal=[1, 0],)
    configs = maxwell_setup([], [s], [m], dx, sz; F, ϵmin)
    @unpack dx, dt, sz, geometry_padding, geometry_staggering, field_padding, source_instances, monitor_instances, u0, = configs

    A = support(monitor_instances[1])
    p = apply(geometry_padding; ϵ, μ, σ, σm)
    p = apply(geometry_staggering, p)
    T1 = l * sqrt(ϵ1) + 2
    Δ = 2
    T2 = T1 + Δ

    u = reduce((u, t) -> maxwell_update!(u, p, t, dx, dt, field_padding, source_instances;), 0:dt:T1, init=deepcopy(u0))
    u, fp, = reduce(T1+dt:dt:T2, init=(u, F(0),)) do (u, fp,), t
        (maxwell_update(u, p, t, dx, dt, field_padding, source_instances),
            fp + dt * flux(u, monitor_instances[1],) / Δ)
    end

    port_powers = mean(fp) .* A
    port_powers, fp
end

function collapse_mode(m, ϵ)
    @unpack Ex, Ey, Ez, Hx, Hy, Hz = m
    Ey, Ex, Hz = map([Ex, Ez, Hy]) do a
        mean(a, dims=2)
    end
    E = sqrt.(abs2.(Ex) + abs2.(Ey) + abs2.(Ez))
    (; Ex, Ey, Hz), sum(E .* ϵ, dims=2) ./ sum(E, dims=2)
end