"""
    function step!(u1, u, p, t, field_padding, source_instances)

Updates fields for 3d. Writes new field into buffer array u1
"""
function step!(u1, u, p, t, dx, dt, field_padding, source_instances)
    ∇ = StaggeredDel([dx, dx, dx])
    ϵ, μ, σ, σm = p
    E, H = u
    E1, H1 = u1
    J = apply(source_instances, t; Jx=zeros(F, size(E[1])), Jy=zeros(F, size(E[2])), Jz=zeros(F, size(E[3])))
    # J = apply(source_instances, t; Jx=0, Jy=0, Jz=0)

    # first update E
    Hx, Hy, Hz = H
    H = mark(field_padding; Hx, Hy, Hz)
    dEdt = (∇ × H - σ * E - J) / ϵ

    E_ = bufferfrom.(E1)
    for i = 1:3
        E_[i][:, :, :] = E[i] + dEdt[i] * dt
    end
    Ex, Ey, Ez = [copy(E_[1]), copy(E_[2]), copy(E_[3])]
    # E .= E + dEdt * dt

    apply!(field_padding; Ex, Ey, Ez)
    H = collect.(H)

    # then update H
    E = mark(field_padding; Ex, Ey, Ez)
    dHdt = -(∇ × E .+ σm * H) / μ

    H_ = bufferfrom.(H1)
    for i = 1:3
        H_[i][:, :, :] = H[i] + dHdt[i] * dt
    end
    Hx, Hy, Hz = [copy(H_[1]), copy(H_[2]), copy(H_[3])]
    # H .= H + dHdt * dt

    apply!(field_padding; Hx, Hy, Hz)
    # E = collect(collect.(E))
    # H = [Hx, Hy, Hz]

    # [E, H]
    u1
end
step3! = step!
# Flux.trainable(m::PaddedArray) = (; a=m.a)

"""
    function step1(u, p, t, field_padding, source_instances)

Updates fields for 1D (Ez, Hy)
"""
function step1(u, p, t, field_padding, source_instances)
    @unpack dx, dt, field_padding, source_instances = configs
    ∇ = Del([dx])
    ϵ, μ, σ, σm = [[a] for a = p]
    Ez, Hy, = u
    Jz, = apply(source_instances, t; Jz=0)

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
    function stepTMz(u, p, t, field_padding, source_instances)

Updates fields for 2d TMz
"""
function stepTMz(u, p, t, field_padding, source_instances)
    @unpack dx, dt, field_padding, source_instances = configs
    ∇ = Del([dx, dx])
    ϵ, μ, σ, σm = [[a] for a = p]
    Ez, Hx, Hy, = u
    Jz, = apply(source_instances, t; Jz=0)

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
    function stepTEz(u, p, t, field_padding, source_instances)

Updates fields for 2d TEz (Hz, Ex, Ey)
"""
function stepTEz(u, p, t, field_padding, source_instances)
    @unpack dx, dt, field_padding, source_instances = configs
    ∇ = Del([dx, dx])
    ϵ, μ, σ, σm = [[a] for a = p]
    Hz, Ex, Ey = u
    Jx, Jy = apply(source_instances, t; Jx=0, Jy=0)

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
