
"""
    function update(u, p, t, field_padding, source_instances; autodiff=true)

Updates fields. mutating if autodiff is false
"""
function update(u, p, t, dx, dt, field_padding, source_instances; autodiff=false, compression=false)
    ϵ, μ, σ, m = group.((p,), (:ϵ, :μ, :σ, :m))
    E, H, J = group.((u,), (:E, :H, :J))
    J0 = J

    J = apply(source_instances, t, J0)
    d = ndims(E(1))
    ∇ = StaggeredDel(fill(dx, d); autodiff)

    # first update E
    H = mark(field_padding, H)
    dEdt = (∇ × H - E * σ - J) / ϵ

    if autodiff
        E = E + dEdt * dt
    else
        for (a, b) = zip(Porcupine.values(E), Porcupine.values(dEdt))
            a .= a + b * dt
        end
    end

    # E = apply_field_padding(field_padding, E; nonzero_only=true)
    H = unmark(H)

    # then update H
    E = mark(field_padding, E)
    dHdt = -(∇ × E + m * H) / μ
    if autodiff
        H = H + dHdt * dt
    else
        for (a, b) = zip(Porcupine.values(H), Porcupine.values(dHdt))
            a .= a + b * dt
        end
    end

    # H = apply(field_padding, H)
    # H = apply_field_padding(field_padding, H; nonzero_only=true)
    # global asdffsd1=H
    E = unmark(E)

    # [E, H]
    if compression
        # return namedtuple([:E => E, :H => H, :J => J0]) 
        # return namedtuple([:E => E, :H => H, :J => J0])
        # return namedtuple([(:E, E), (:H, H), (:J, J0)])
        # return namedtuple([[:E, E], [:H, H], [:J, J0]])
        return (E, H, J)
    end

    # (E, H, J0)
    merge(E, H, J0)
    # merge(E + σ, H + μ, J0)
end
