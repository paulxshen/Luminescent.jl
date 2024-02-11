"""
    function step3(u, p, t, field_padding, source_effects)

Updates fields for 3d
"""
function step3!(u, p, t, dx, dt, field_padding, source_effects)
    ∇ = StaggeredDel([dx, dx, dx])
    ϵ, μ, σ, σm = p
    E, H = u
    J = apply(source_effects, t; Jx=0, Jy=0, Jz=0)

    # first update E
    Hx, Hy, Hz = H
    H = mark(field_padding; Hx, Hy, Hz)
    dEdt = (∇ × H - σ * E - J) / ϵ

    E_ = bufferfrom.(E)
    for i = 1:3
        E_[i][:, :, :] = E[i] + dEdt[i] * dt
    end
    Ex, Ey, Ez = [copy(E_[1]), copy(E_[2]), copy(E_[3])]
    # E .= E + dEdt * dt

    Ex, Ey, Ez = E
    apply!(field_padding; Ex, Ey, Ez)
    H = collect.(H)

    # then update H
    E = mark(field_padding; Ex, Ey, Ez)
    dHdt = -(∇ × E .+ σm * H) / μ

    H_ = bufferfrom.(H)
    for i = 1:3
        H_[i][:, :, :] = H[i] + dHdt[i] * dt
    end
    Hx, Hy, Hz = [copy(H_[1]), copy(H_[2]), copy(H_[3])]
    # H .= H + dHdt * dt

    Hx, Hy, Hz = H
    apply!(field_padding; Hx, Hy, Hz)
    E = collect.(E)
    H = Hx, Hy, Hz

    [E, H]
end
step! = step3!
# Flux.trainable(m::PaddedArray) = (; a=m.a)

"""
    function step1(u, p, t, field_padding, source_effects)

Updates fields for 1D (Ez, Hy)
"""
function step1(u, p, t, field_padding, source_effects)
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
    function stepTMz(u, p, t, field_padding, source_effects)

Updates fields for 2d TMz
"""
function stepTMz(u, p, t, field_padding, source_effects)
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
    function stepTEz(u, p, t, field_padding, source_effects)

Updates fields for 2d TEz (Hz, Ex, Ey)
"""
function stepTEz(u, p, t, field_padding, source_effects)
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
