# Base.-()
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
"""
    function step3(u, p, t, field_padding, source_effects)

Updates fields for 3d
"""
function step3!(u, p, t, dx, dt, field_padding, source_effects; Buffer=nothing, bufferfrom=nothing)
    ∇ = StaggeredDel([dx, dx, dx])
    ϵ, μ, σ, σm = p
    Ex, Ey, Ez, Hx, Hy, Hz = u
    J = apply(source_effects, t; Jx=0, Jy=0, Jz=0)
    autodiff = !isnothing(bufferfrom)

    # first update E
    H = mark(field_padding; Hx, Hy, Hz)
    E = [Ex, Ey, Ez]

    dEdt = (∇ × H .- σ * E .- J) / ϵ
    # println(typeof.((J, ϵ)))
    H = collect.(H)
    E_ = autodiff ? bufferfrom.(E) : E
    for i = 1:3
        # if autodiff
        #     E_ = Buffer.(E)
        # end
        E_[i][:, :, :] = E[i] + dEdt[i] * dt
        # E_[i] .= E[i] + dEdt[i] * dt
        # E = E + dEdt * dt
    end

    Ex, Ey, Ez = E_
    apply(field_padding; Ex, Ey, Ez)
    E = autodiff ? (copy.(E_)) : E_
    # println(extrema(E[3]))

    # then update H
    Ex, Ey, Ez = E
    E = mark(field_padding; Ex, Ey, Ez)

    dHdt = -(∇ × E .+ σm * H) / μ
    E = collect.(E)
    H_ = autodiff ? bufferfrom.(H) : H
    for i = 1:3
        # if autodiff
        #     H_ = Buffer.(H)
        # end
        H_[i][:, :, :] = H[i] + dHdt[i] * dt
        # H_[i] .= H[i] + dHdt[i] * dt
        # H_ = H + dHdt * dt
    end

    #     H = copy.(H)
    Hx, Hy, Hz = H_
    apply(field_padding; Hx, Hy, Hz)
    # A = autodiff ? (copy.(H_)) : H_

    # H = collect.(H)
    if autodiff
        return [E..., copy(Hx), copy(Hy), copy(Hz)]
    end
    [E..., H...]
    # [E[1], E[2], E[3], H[1], H[2], H[3]]
end
step! = step3!
# Flux.trainable(m::PaddedArray) = (; a=m.a)