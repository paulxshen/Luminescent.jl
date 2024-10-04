
"""
    function update(u, p, t, field_padding, source_instances; autodiff=true)

Updates fields. mutating if autodiff is false
"""
function update(u, p, t, dx, dt, field_padding, source_instances; past=false, alg=nothing)
    # t, dx, dt, field_padding, source_instances, autodiff, save_memory = ignore_derivatives() do
    #     t, dx, dt, field_padding, source_instances, autodiff, save_memory
    # end
    if past
        _u, _dudt, u = u
    end

    ϵ, μ, σ, m = group.((p,), (:ϵ, :μ, :σ, :m))
    E, H, J = group.((u,), (:E, :H, :J))
    J0 = J
    T = eltype(t)

    J = apply(source_instances, t, J0)
    d = ndims(E(1))
    padamt = ignore_derivatives() do
        onedge = namedtuple([k => [sum(left, v) sum(right, v)] for (k, v) = pairs(field_padding)])
        1 - onedge
    end
    ∇ = StaggeredDel(fill(dx, d), padamt, alg)

    # first update E
    dEdt = (∇ × H - E * σ - J) / ϵ

    if !past
        E = E + dEdt * dt
    else
        # dt=1
        ts = dt * range(-2, -0.5, 4)
        f(t) = t .^ (0:3)
        f_(t) = [0, 1, 2t, 3t^2]
        C = inv(hcat(f(ts[1]), f_(ts[2]), f(ts[3]), f_(ts[4]))')[1, :]
        C = T.(C)

        _E, _H = group.((_u,), (:E, :H))
        _dEdt, _dHdt = _dudt

        _E = deepcopy(E)
        E = sum(C .* [_E, _dEdt, E, dEdt])
    end
    # else
    #     for (a, b) = zip(Porcupine.values(E), Porcupine.values(dEdt))
    #         a .= a + b * dt
    #     end
    # end

    # E = apply_field_padding(field_padding, E; nonzero_only=true)

    # then update H
    dHdt = -(∇ × E + m * H) / μ
    if !past
        H = H + dHdt * dt
    else
        _H = deepcopy(H)
        H = sum(C .* [_H, _dHdt, H, dHdt])
    end
    #     for (a, b) = zip(Porcupine.values(H), Porcupine.values(dHdt))
    # a .= a + b * dt
    # end

    # H = apply(field_padding, H)
    # H = apply_field_padding(field_padding, H; nonzero_only=true)
    # global asdffsd1=H

    u = merge(E, H, J0)
    if !past
        return u
    end
    _dudt = (dEdt, dHdt)
    (_u, _dudt, u)
end
