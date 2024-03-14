
"""
    function maxwell_update!(u, p, t, field_padding, source_instances)
    function maxwell_update!(u1, u, p, t, field_padding, source_instances)

Updates fields for 3d. Please use `maxwell_update` instead of `maxwell_update!` when doing AD. Mutating `maxwell_update!` Writes new fields either onto old fields or into buffer arrays u1
"""
function maxwell_update!(u1, u, p, t, dx, dt, field_padding, source_instances)
    ∇ = StaggeredDel([dx, dx, dx])
    @unpack ϵ, μ, σ, σm = p
    # ϵ, μ, σ, σm = p |> values
    @unpack E, H = u
    E1 = u1[:E]
    H1 = u1[:H]
    # J = apply(source_instances, t; Jx=zeros(F, size(E[1])), Jy=zeros(F, size(E[2])), Jz=zeros(F, size(E[3])))
    J = apply(source_instances, t; Jx=0, Jy=0, Jz=0)

    # first update E
    Hx, Hy, Hz = H
    H = mark(field_padding; Hx, Hy, Hz)
    dEdt = (∇ × H - σ * E - J) / ϵ

    for i = 1:3
        E1[i][:, :, :] = E[i] + dEdt[i] * dt
    end
    Ex, Ey, Ez = E1
    apply!(field_padding; Ex, Ey, Ez)
    H = unmark(; H...)

    # then update H
    E = mark(field_padding; Ex, Ey, Ez)
    dHdt = -(∇ × E .+ σm * H) / μ

    for i = 1:3
        H1[i][:, :, :] = H[i] + dHdt[i] * dt
    end
    Hx, Hy, Hz = H1
    apply!(field_padding; Hx, Hy, Hz)
    E = unmark(; E...)

    [E, H]
    # u1
end
maxwell_update!(u, p, t, dx, dt, field_padding, source_instances) = maxwell_update!(u, u, p, t, dx, dt, field_padding, source_instances)

"""
    function maxwell_update(u, p, t, field_padding, source_instances)

Updates fields for 3d in a manner amenable to AD. See also Mutating `maxwell_update!`
"""
function maxwell_update(u, p, t, dx, dt, field_padding, source_instances; ignore_boundary_autodiff=false)
    ∇ = StaggeredDel([dx, dx, dx])
    @unpack ϵ, μ, σ, σm = p
    @unpack E, H = u
    J = apply(source_instances, t; Jx=0, Jy=0, Jz=0)

    # first update E
    # H = mark(field_padding; H...)
    H = mark(field_padding, H)
    # @info typeof(E)
    dEdt = (∇ × H - σ * E - J) / ϵ
    # @info typeof(dEdt)
    E = E + dEdt * dt

    # Ex, Ey, Ez = E
    # E = apply(field_padding; Ex, Ey, Ez)
    # Hx, Hy, Hz = H
    # H = unmark(; Hx, Hy, Hz)

    # E = apply(field_padding; E...)
    # H = unmark(; H...)
    E = apply(field_padding, E)
    H = unmark(H)

    # then update H
    # E = mark(field_padding; Ex, Ey, Ez)
    E = mark(field_padding, E)
    dHdt = -(∇ × E .+ σm * H) / μ
    H = H + dHdt * dt

    # Ex, Ey, Ez = E
    # Ex, Ey, Ez = unmark(; Ex, Ey, Ez)
    # E = (; Ex, Ey, Ez)
    # Hx, Hy, Hz = H
    # H = apply(field_padding; Hx, Hy, Hz)

    # H = apply(field_padding; H...)
    # E = unmark(; E...)
    H = apply(field_padding, H)
    E = unmark(E)

    # [E, H]
    MyDict([:E => E, :H => H])
    # (; E, H)
end
maxwell_update! = maxwell_update!
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
