function calibrate_mode(mode, ϵmode, dx; F=Float32, verbose=false, name="")
    mode, ϵmode, dx = (mode, ϵmode, dx) |> F
    d = ndims(ϵmode) + 1
    n = sqrt(mean(ϵmode))
    l = 1.6
    Δ = F[l*n+2, 2]
    T = cumsum(Δ)

    # wm = wwg / 2
    # w = 2wm + wwg
    # sz = (round(Int, l / dx), round(Int, w / dx))
    ϵmin = minimum(ϵmode)
    L = size(ϵmode) * dx
    _margin = round(0.5 / dx)
    margin = _margin * dx
    Lm += 2margin
    sz = round.(Int, (l / dx, (Lm / dx)...))

    μ = ones(F, sz)
    σ = zeros(F, sz)
    σm = zeros(F, sz)

    ϵ = stack(fill(ϵmode, sz[1]))
    ϵ = pad(ϵ, :replicate, (0, _margin), (0, _margin))
    # ϵ |> heatmap |> display
    # ϵ[1+round(wm / dx):end-round(wm / dx)] .= ϵ2

    normal = [1, zeros(d - 1)]
    m = [ModalMonitor(mode, [x, Lm / 2...], L, normal,) for x = 0.5:0.5:1.5]
    s = ModalSource(t -> cispi(2t), mode, [margin_offset, Lm / 2...], L, normal;)
    prob = maxwell_setup([], [s], m, dx, sz; F, ϵmin, verbose)
    # recordsim("$(@__DIR__)/$name.mp4", h, ;
    #     dt,
    #     field=:Hz,
    #     monitor_instances,
    #     source_instances,
    #     geometry=p.ϵ.ϵyy,
    # )

    @show port_powers = mean.(pp) .* A
    fp[end], a[end][1], pp
end
