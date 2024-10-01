
"""
    function update(u, p, t, field_padding, source_instances; autodiff=true)

Updates fields. mutating if autodiff is false
"""
function update(u, p, t, dx, dt, field_padding, source_instances; autodiff=false)
    # t, dx, dt, field_padding, source_instances, autodiff, save_memory = ignore_derivatives() do
    #     t, dx, dt, field_padding, source_instances, autodiff, save_memory
    # end
    ϵ, μ, σ, m = group.((p,), (:ϵ, :μ, :σ, :m))
    E, H, J = group.((u,), (:E, :H, :J))
    J0 = J

    J = apply(source_instances, t, J0)
    d = ndims(E(1))
    padamt = ignore_derivatives() do
        onedge = namedtuple([k => [sum(left, v) sum(right, v)] for (k, v) = pairs(field_padding)])
        1 - onedge
    end
    ∇ = StaggeredDel(fill(dx, d), padamt)

    # first update E
    dEdt = (∇ × H - E * σ - J) / ϵ

    # if autodiff
    E = E + dEdt * dt
    # else
    #     for (a, b) = zip(Porcupine.values(E), Porcupine.values(dEdt))
    #         a .= a + b * dt
    #     end
    # end

    # E = apply_field_padding(field_padding, E; nonzero_only=true)

    # then update H
    dHdt = -(∇ × E + m * H) / μ
    # if autodiff
    H = H + dHdt * dt
    # else
    #     for (a, b) = zip(Porcupine.values(H), Porcupine.values(dHdt))
    # a .= a + b * dt
    #     end
    # end

    # H = apply(field_padding, H)
    # H = apply_field_padding(field_padding, H; nonzero_only=true)
    # global asdffsd1=H

    merge(E, H, J0)
    # merge(E + σ, H + μ, J0)
end
