
"""
    function update(u, p, t, field_padding, source_instances; autodiff=true)

Updates fields. mutating if autodiff is false
"""
function update(u, p, t, dx, dt, field_padding, source_instances; autodiff=false, save_memory=false)
    @unpack ϵ, μ, σ, m = p
    # ϵ, μ, σ, m = Ref.([ϵ, μ, σ, m])
    E, H, J0 = unroll(u, :E, :H, :J)
    # J0 = deepcopy(J)
    J = apply(source_instances, t, J0)
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

    E = apply(field_padding, E; nonzero_only=true)
    # E = apply(field_padding, E)
    H = unmark(H)

    # then update H
    E = mark(field_padding, E)
    dHdt = -(∇ × E + m * H) / μ
    if autodiff
        H = H + dHdt * dt
    else
        for (a, b) = zip(values(H), values(dHdt))
            a .= a + b * dt
        end
    end

    # H = apply(field_padding, H)
    H = apply(field_padding, H; nonzero_only=true)
    E = unmark(E)

    # [E, H]
    if save_memory
        # return NamedTuple([:E => E, :H => H, :J => J0]) 
        # return namedtuple([:E => E, :H => H, :J => J0])
        # return namedtuple([(:E, E), (:H, H), (:J, J0)])
        # return namedtuple([[:E, E], [:H, H], [:J, J0]])
        return (E, H, J)
    end
    (; E, H, J=J0)
end
