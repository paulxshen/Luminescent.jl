
"""
    function update(u, p, t, field_padding, source_instances; autodiff=true)

Updates fields. mutating if autodiff is false
"""
function update(u, p, t, dx, dt, field_padding, source_instances; autodiff=true)
    @unpack ϵ, μ, σ, σm = p
    # ϵ, μ, σ, σm = Ref.([ϵ, μ, σ, σm])
    @unpack E, H, J = u
    J0 = deepcopy(J)
    J = apply(source_instances, t, J)
    d = ndims(first(values(E)))
    ∇ = StaggeredDel(fill(dx, d); autodiff)

    # first update E
    H = mark(field_padding, H)
    dEdt = (∇ × H - σ * E - J) / ϵ

    if autodiff
        E = E + dEdt * dt
    else
        for (a, b) = zip(values(E), values(dEdt))
            a .= a + b * dt
        end
    end

    E = apply(field_padding, E)
    H = unmark(H)

    # then update H
    E = mark(field_padding, E)
    dHdt = -(∇ × E + σm * H) / μ
    if autodiff
        H = H + dHdt * dt
    else
        for (a, b) = zip(values(H), values(dHdt))
            a .= a + b * dt
        end
    end

    H = apply(field_padding, H)
    E = unmark(E)

    # [E, H]
    # dict([:E => E, :H => H, :J => J0])
    (; E, H, J=J0)
end
