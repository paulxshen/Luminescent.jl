function calibrate_mode(mode, ϵmode, dx, ; F=Float32, verbose=false, name="", kwargs...)
    mode, ϵmode, dx, = (mode, ϵmode, dx,) |> F
    d = ndims(ϵmode) + 1

    source_margin, port_source_offset, mode_margin = trim.((SOURCE_MARGIN, PORT_SOURCE_OFFSET, MODE_MARGIN), dx)
    l0 = 2
    l = l0 + port_source_offset + source_margin

    # wm = wwg / 2
    # w = 2wm + wwg
    # sz = (round(Int, l / dx), round(Int, w / dx))
    L = [size(ϵmode) * dx...]
    sz = (round(l / dx), size(ϵmode)...)

    ϵ = stack(fill(ϵmode, sz[1]))
    if d == 2
        ϵ = ϵ'
        ϵ = pad(ϵ, :replicate, (0, round(mode_margin / dx)))
    else
        ϵ = permutedims(ϵ, (3, 1, 2))
    end
    sz = size(ϵ)
    l, w, = sz * dx

    normal = [1, zeros(d - 1)...]
    tangent = [0, -1, zeros(d - 2)...]
    m = [Monitor(mode, [x, w / 2], normal, tangent, L) for x = range(0, l0, 3) + source_margin + port_source_offset]
    s = Source(t -> cispi(2t), mode, [source_margin, w / 2], -normal, tangent, L;)

    prob = setup([], [s], m, dx, sz; F, ϵ, verbose)
    sol = solve(prob)

    @unpack fields, modes, forward_mode_powers = sol
    mode = modes[1][1]
    power = forward_mode_powers[1][1][1]
    @show forward_mode_powers
    quickie(sol) |> display
    global ret = (; mode, power, sol)

end
