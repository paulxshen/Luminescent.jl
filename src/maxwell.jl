"""
    function step1(u, p, t, configs)

Updates fields for 1D (Ez, Hy)
"""
function step1(u, p, t, configs)
    @unpack dx, dt, field_padding, source_effects = configs
    ∇ = Del([dx])
    ϵ, μ, σ, σm = [[a] for a = p]
    Ez, Hy, = u
    Jz, = apply(source_effects, t; Jz=0)

    # first update E
    H_ = apply(field_padding; Hy=PaddedArray(Hy))

    E = [Ez]
    J = [Jz]

    dEdt = (∇ × H_ - σ * E - J) / ϵ
    E += dEdt * dt
    Ez, = E

    # then update H
    Ez_, = apply(field_padding; Ez=PaddedArray(Ez))
    E_ = [Ez_]
    H = [Hy]
    dHdt = -(∇ × E_ + σm * H) / μ
    H += dHdt * dt
    Hy, = H

    [Ez, Hy,]
end
"""
    function stepTMz(u, p, t, configs)

Updates fields for 2d TMz
"""
function stepTMz(u, p, t, configs)
    @unpack dx, dt, field_padding, source_effects = configs
    ∇ = Del([dx, dx])
    ϵ, μ, σ, σm = [[a] for a = p]
    Ez, Hx, Hy, = u
    Jz, = apply(source_effects, t; Jz=0)

    # first update E
    Hx_, Hy_ = apply(field_padding; Hx=PaddedArray(Hx), Hy=PaddedArray(Hy))
    H_ = [Hx_, Hy_]

    E = [Ez]
    J = [Jz]

    dEdt = (∇ × H_ - σ * E - J) / ϵ
    E += dEdt * dt
    Ez, = E

    # then update H
    E_ = apply(field_padding; Ez=PaddedArray(Ez))
    H = [Hx, Hy]
    dHdt = -(∇ × E_ + σm * H) / μ
    H += dHdt * dt
    Hx, Hy = H

    [Ez, Hx, Hy,]
end
"""
    function stepTEz(u, p, t, configs)

Updates fields for 2d TEz (Hz, Ex, Ey)
"""
function stepTEz(u, p, t, configs)
    @unpack dx, dt, field_padding, source_effects = configs
    ∇ = Del([dx, dx])
    ϵ, μ, σ, σm = [[a] for a = p]
    Hz, Ex, Ey = u
    Jx, Jy = apply(source_effects, t; Jx=0, Jy=0)

    # first update E
    H_ = apply(field_padding; Hz=PaddedArray(Hz))


    E = [Ex, Ey]
    J = [Jx, Jy]

    dEdt = (∇ × H_ - σ * E - J) / ϵ
    E += dEdt * dt
    Ex, Ey = E

    # then update H
    E_ = apply(field_padding; Ex=PaddedArray(Ex), Ey=PaddedArray(Ey))
    H = [Hz]
    dHdt = -(∇ × E_ + σm * H) / μ
    H += dHdt * dt
    Hz, = H

    [Hz, Ex, Ey]
end
"""
    function step3(u, p, t, configs)

Updates fields for 3d
"""
function step3(u, p, t, configs)
    @unpack dx, dt, field_padding, source_effects = configs
    ∇ = Del([dx, dx, dx])
    ϵ, μ, σ, σm = [[a] for a = p]
    Ex, Ey, Ez, Hx, Hy, Hz = u
    Jx, Jy, Jz, = apply(source_effects, t; Jx=0, Jy=0, Jz=0)

    # first update E
    Hx_, Hy_, Hz_ = apply(field_padding; Hx=PaddedArray(Hx), Hy=PaddedArray(Hy), Hz=PaddedArray(Hz))
    H_ = [Hx_, Hy_, Hz_]

    E = [Ex, Ey, Ez]
    J = [Jx, Jy, Jz]

    dEdt = (∇ × H_ - σ * E .- J) / ϵ
    E += dEdt * dt
    Ex, Ey, Ez, = E

    # then update H
    E_ = apply(field_padding; Ex=PaddedArray(Ex), Ey=PaddedArray(Ey), Ez=PaddedArray(Ez))
    H = [Hx, Hy, Hz]
    dHdt = -(∇ × E_ + σm * H) / μ
    H += dHdt * dt
    Hx, Hy, Hz = H

    [Ex, Ey, Ez, Hx, Hy, Hz]
end

# Flux.trainable(m::PaddedArray) = (; a=m.a)