"""

Updates fields. 
"""
function update(u, p, t, dt, field_diffdeltas, field_diffpadvals, source_instances; alg=nothing)
    # unpack fields and geometry
    @unpack E, H = u
    @unpack μ, m, invϵ = p

    # dispersive material parameters
    poles = @ignore_derivatives sort(Set([k[2:end-1] for k = string.(keys(u)) if startswith(k, "J")]))

    Jkeys = Symbol.(("J",) .* poles)
    Jm = [u(k) for k = Jkeys]
    Pkeys = Symbol.(("P",) .* poles)
    Pm = [u(k) for k = Pkeys]

    σm = [p("σ$k") for k = poles]
    γm = [p("γ$k") for k = poles]
    βm = [p("β$k") for k = poles]

    @nogradvars t, dt, field_diffdeltas, field_diffpadvals, source_instances, μ, m
    N = ndims(E(1))

    # staggered grid housekeeping
    ∇ = Del(field_diffdeltas, field_diffpadvals)

    # inject sources
    Js = sum([s(t) for s = source_instances])

    # update time dependent polarization and displacement fields for dispersive materials
    Jm = [
        if all(==(0), γ)
            J = 0J + E ⊙ σ
        else
            dJdt = γ ⊙ (E ⊙ σ - J - β ⊙ P)
            J += dJdt * dt
        end for (J, σ, γ, β) = zip(Jm, σm, γm, βm)
    ]
    Pm = [
        if any(!=(0), β)
            P += J * dt
        else
            P
        end for (P, J, β) = zip(Pm, Jm, βm)
    ]


    # first update E
    # tensor subpixel smoothing
    # dEdt = (∇ × H - Js - sum(Jm)) ⊘ ϵ
    dEdt = invϵ * (∇ × H - Js - sum(Jm))
    E += dEdt * dt

    # then update H
    dHdt = -(∇ × E + H ⊙ m) ⊘ μ
    H += dHdt * dt

    @ignore_derivatives gc()
    # @ignore_derivatives unsafe_free!.((Js, dEdt, dHdt))
    # (; E, H, (Jkeys .=> Jm)..., (Pkeys .=> Pm)...)
    # namedtuple([:E => E, :H => H, (Jkeys .=> Jm)..., (Pkeys .=> Pm)...])
    namedtuple(vcat([:E => E, :H => H], Jkeys .=> Jm, Pkeys .=> Pm))
end
