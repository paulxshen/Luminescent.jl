function step_TMz(u, p, t, configs)
    @unpack dx, dt, field_padding, source_effects = configs
    ∇ = Del([dx, dx])
    ϵ, μ, σ, σm = [[a] for a = p]
    Ez, Hx, Hy, = u
    Jz, = apply(source_effects, t; Jz=0)

    # first update E
    H_ = apply(field_padding, ; Hx=PaddedArray(Hx), Hy=PaddedArray(Hy))
    E = [Ez]
    J = [Jz]
    dEdt = (∇ × H_ - σ * E - J) / ϵ
    E += dEdt * dt
    Ez, = E

    # then update H
    E_ = apply(field_padding, ; Ez=PaddedArray(Ez))
    H = [Hx, Hy]
    dHdt = -(∇ × E_ + σm * H) / μ
    H += dHdt * dt
    Hx, Hy = H

    [Ez, Hx, Hy,]
end