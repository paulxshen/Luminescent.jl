
"""
    function update!(u, p, t, field_padding, source_instances)
    function update!(u1, u, p, t, field_padding, source_instances)

Updates fields for 3d. Please use `update` instead of `update!` when doing AD. Mutating `update!` Writes new fields either onto old fields or into buffer arrays u1
"""
function update!(u1, u, p, t, dx, dt, field_padding, source_instances)
    @unpack ϵ, μ, σ, σm = p
    # ϵ, μ, σ, σm = Ref.([ϵ, μ, σ, σm])

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
update!(u, p, t, dx, dt, field_padding, source_instances) = update!(u, u, p, t, dx, dt, field_padding, source_instances)

"""
    function update(u, p, t, field_padding, source_instances)
    
Updates fields for 3d in a manner amenable to AD. See also Mutating `update!`
        """
function update(u, p, t, dx, dt, field_padding, source_instances)
    @unpack ϵ, μ, σ, σm = p
    # ϵ, μ, σ, σm = Ref.([ϵ, μ, σ, σm])
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
    # @show extrema((values(dEdt))[2]), extrema(J.Jy)
    E = E + dEdt * dt

    E = apply(field_padding, E)
    H = unmark(H)

    # then update H
    E = mark(field_padding, E)
    dHdt = -(∇ × E + σm * H) / μ
    H = H + dHdt * dt

    H = apply(field_padding, H)
    E = unmark(E)

    # [E, H]
    # dict([:E => E, :H => H, :J => J0])
    (; E, H, J=J0)
end
update! = update!
# Flux.trainable(m::PaddedArray) = (; a=m.a)

