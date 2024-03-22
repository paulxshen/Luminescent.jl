
"""
    function maxwell_update!(u, p, t, field_padding, source_instances)
    function maxwell_update!(u1, u, p, t, field_padding, source_instances)

Updates fields for 3d. Please use `maxwell_update` instead of `maxwell_update!` when doing AD. Mutating `maxwell_update!` Writes new fields either onto old fields or into buffer arrays u1
"""
function maxwell_update!(u1, u, p, t, dx, dt, field_padding, source_instances)
    @unpack ϵ, μ, σ, σm = p
    @unpack E, H, J = u
    J0 = deepcopy(J)
    E1 = u1[:E]
    H1 = u1[:H]
    J = apply(source_instances, t; J...)
    d = ndims(first(values(E)))
    ∇ = StaggeredDel(fill(dx, d))

    # first update E
    H = mark(field_padding, H)
    dEdt = (∇ × H - σ * E - J) / ϵ

    for (a, b, c) = zip(values(E1), values(E), values(dEdt))
        a .= b + c * dt
    end
    E = apply!(field_padding, E1)
    H = unmark(H)

    # then update H
    E = mark(field_padding, E)
    dHdt = -(∇ × E + σm * H) / μ

    for (a, b, c) = zip(values(H1), values(H), values(dHdt))
        a .= b + c * dt
    end
    H = apply!(field_padding, H1)
    E = unmark(E)

    (; E, H, J=J0)
    # u1
end
maxwell_update!(u, p, t, dx, dt, field_padding, source_instances) = maxwell_update!(u, u, p, t, dx, dt, field_padding, source_instances)

"""
    function maxwell_update(u, p, t, field_padding, source_instances)
    
Updates fields for 3d in a manner amenable to AD. See also Mutating `maxwell_update!`
        """
function maxwell_update(u, p, t, dx, dt, field_padding, source_instances; ignore_boundary_autodiff=false)
    @unpack ϵ, μ, σ, σm = p
    @unpack E, H, J = u
    J0 = deepcopy(J)
    J = apply(source_instances, t, J)
    d = ndims(first(values(E)))
    ∇ = StaggeredDel(fill(dx, d))

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
    dHdt = -(∇ × E + σm * H) / μ
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
    # dict([:E => E, :H => H, :J => J0])
    (; E, H, J=J0)
end
maxwell_update! = maxwell_update!
# Flux.trainable(m::PaddedArray) = (; a=m.a)

