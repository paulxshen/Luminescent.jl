"""
    function update(u, p, t, field_padding, source_instances)

Updates fields. 
"""
function update(u, p, t, dx, dt, field_padding, source_instances; alg=nothing)
    # unpack fields and geometry
    E, H, J = [u(Regex("$k.*")) for k = (:E, :H, :J)]
    @unpack invϵ, ϵ, μ, σ, m = p
    J0 = J

    # staggered grid housekeeping
    N = ndims(E(1))
    padamt = ignore_derivatives() do
        1 - namedtuple([
            k => [sum(left, v) sum(right, v)]
            for (k, v) = pairs(field_padding)
        ]) |> cpu
    end
    Epadamt = padamt(r"E.*")
    ∇ = StaggeredDel(fill(dx, N), padamt, alg)

    # inject sources
    J = apply(source_instances, t, J0)

    # first update E
    # dEdt = (∇ × H - E * σ - J) / ϵ
    # dEdt = invϵ * (∇ × H - E * σ - J)
    dDdt = (∇ × H - E * σ - J)
    # tensor subpixel smoothing with staggered grid
    dEdt = map(1:size(invϵ, 1)) do i
        row = invϵ[i, :]
        sum(map(eachindex(row)) do j
            l, r = eachcol((Epadamt[i] - Epadamt[j]) / 2)
            crop(row[j] .* dDdt[j], l, r)
        end)
    end
    E += dEdt * dt

    # then update H
    dHdt = -(∇ × E + m * H) / μ
    H += dHdt * dt

    u = merge(E, H, J0)
end

# E = apply_field_padding(field_padding, E; nonzero_only=true)
# """
#     function update(u, p, t, field_padding, source_instances; autodiff=true)

# Updates fields. mutating if autodiff is false
# """
# function update(u, p, t, dx, dt, field_padding, source_instances; past=false, alg=nothing)
#     # t, dx, dt, field_padding, source_instances, autodiff, save_memory = ignore_derivatives() do
#     #     t, dx, dt, field_padding, source_instances, autodiff, save_memory
#     # end
#     # _u, _dudt, u = u
#     @unpack invϵ, ϵ, μ, σ, m = p
#     E, H, J = group.((u,), (:E, :H, :J))
#     J0 = J
#     Hkeys = keys(H)
#     Ekeys = keys(E)
#     T = eltype(t)

#     J = apply(source_instances, t, J0)
#     d = ndims(E(1))
#     padamt = ignore_derivatives() do
#         onedge = namedtuple([k => [sum(left, v) sum(right, v)] for (k, v) = pairs(field_padding)]) |> cpu
#         1 - onedge
#     end
#     Epadamt = padamt(r"E.*")
#     ∇ = StaggeredDel(fill(dx, d), padamt, alg)

#     # first update E
#     # dEdt = (∇ × H - E * σ - J) / ϵ
#     # dEdt = invϵ * (∇ × H - E * σ - J)
#     dDdt = (∇ × H - E * σ - J)
#     # dEdt = map(Epadamt, eachrow(invϵ)) do Epadamti, row
#     dEdt = map(1:size(invϵ, 1)) do i
#         Epadamti = Epadamt[i]
#         row = invϵ[i, :]
#         # sum(map(Epadamt, row, dDdt) do Epadamtj, a, d
#         sum(map(eachindex(row)) do j
#             Epadamtj = Epadamt[j]
#             a = row[j]
#             d = dDdt[j]
#             l, r = eachcol((Epadamti - Epadamtj) / 2)
#             crop(a .* d, l, r)
#             # a .* d
#         end)
#     end
#     # dEdt = dDdt / 3

#     # dt=.5
#     # C = ignore_derivatives() do

#     #     ts = dt * range(-2, -0.5, 4)
#     #     f(t) = t .^ (0:3)
#     #     f_(t) = [0, 1, 2t, 3t^2]
#     #     C = inv(hcat(f(ts[1]), f_(ts[2]), f(ts[3]), f_(ts[4]))')[1, :]
#     # end
#     # C = T.(C)

#     # _E, _H = group.((_u,), (:E, :H))
#     # _dEdt, _dHdt = group.((_dudt,), (:E, :H))

#     # E1 = deepcopy(E)
#     # E = sum(C .* [_E, _dEdt, E, dEdt])
#     E = E + dEdt * dt

#     # dEdt = namedtuple(Pair.(Ekeys, values(dEdt)))
#     # E = namedtuple(Pair.(Ekeys, values(E)))
#     # else
#     #     for (a, b) = zip(Porcupine.values(E), Porcupine.values(dEdt))
#     #         a .= a + b * dt
#     #     end
#     # end

#     # E = apply_field_padding(field_padding, E; nonzero_only=true)

#     # then update H
#     dHdt = -(∇ × E + m * H) / μ
#     # H1 = deepcopy(H)
#     # H = sum(C .* [_H, _dHdt, H, dHdt])
#     H = H + dHdt * dt

#     # dHdt = namedtuple(Pair.(Hkeys, values(dHdt)))
#     # H = namedtuple(Pair.(Hkeys, values(H)))
#     #     for (a, b) = zip(Porcupine.values(H), Porcupine.values(dHdt))
#     # a .= a + b * dt
#     # end

#     # H = apply(field_padding, H)
#     # H = apply_field_padding(field_padding, H; nonzero_only=true)
#     # global asdffsd1=H

#     u = merge(E, H, J0)
#     # _u = merge(E1, H1)
#     # _dudt = merge(dEdt, dHdt)
#     # (_u, _dudt, u)
# end
