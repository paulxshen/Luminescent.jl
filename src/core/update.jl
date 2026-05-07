function update(u, p, t; dt, ∇ₕ, ∇ₑ, source_instances)
    # unpack fields and geometry
    # global _p = p
    @unpack E, H, Pn, Jn = u
    @unpack m, σ, invϵ, invμ, Cn, χ2, χ3 = p
    @nograd t, dt, ∇ₕ, ∇ₑ, source_instances, invμ, m, Cn, χ2, χ3
    @nograd σ, Pn, Jn #nc

    N = length(E)
    # global _u, _p = u, p
    # error("debug update")
    # inject sources
    Js = @ignore_derivatives sum(source_instances) do s
        s(t)
    end

    if Cn != 0
        # https://meep.readthedocs.io/en/latest/Materials/#material-dispersion
        Jn += map(Pn, Jn, Cn) do P, J, (γ, ω2, islorentzian)
            (ω2 ⊙ (E ⊙ σ - islorentzian ⊙ P) - γ ⊙ J) * dt
        end
        Pn += dt * Jn
    end

    if χ2 != 0 || χ3 != 0
        dϵ = χ2 ⊙ E + χ3 ⊙ [Porcupine.norm2(E)]
        invϵ = map(Iterators.product(1:N, 1:N)) do (i, j)
            a = invϵ[i, j]
            i == j ? 1 ⊘ (1 ⊘ a + dϵ[i]) : a
        end
    end

    # first update E
    # tensor subpixel smoothing
    # dDdt = ∇ₕ × H - Js - sum(Jn) - σ ⊙ E
    # dDdt = ∇ₕ × H - Js - σ ⊙ E
    # _H = 𝐟(H)
    # _E = 𝐟(map(_values(E)))
    # dDdt = (∇ₕ × _H - Js - map(_E) do a
    #     a .* σ[1] #nc
    # end)
    dDdt = (∇ₕ × H - Js - map(_values(E)) do a
        a .* σ[1] #nc
    end)
    dEdt = invϵ * dDdt
    E += dEdt * dt

    # then update H
    # dBdt = -(∇ₑ × E + H ⊙ m)
    # _E = 𝐟(E)
    # _H = 𝐟(map(_values(H)))
    # dBdt = -(∇ₑ × _E + map(_H) do a
    #     a .* m[1]
    # end)
    dBdt = -(∇ₑ × E + map(_values(H)) do a
        a .* m[1]
    end)
    dHdt = invμ ⊙ dBdt
    H += dHdt * dt
    # @ignore_derivatives GC.gc()

    (; E, H, Pn, Jn)
end
