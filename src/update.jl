"""
    function update(u, p, t, field_boundvals, source_instances)

Updates fields. 
"""
function update(u, p, t, dt, field_diffdeltas, field_diffpadvals, source_instances; alg=nothing)
    # unpack fields and geometry
    E, H = [u(ignore_derivatives() do
        Regex("$k.*")
    end) for k = (:E, :H, :J)]
    # global ___p = p
    @unpack invϵ, μ, σ, m, ϵ = p
    T = eltype(t)

    # staggered grid housekeeping
    N = ndims(E(1))
    ∇ = ignore_derivatives() do
        Del(field_diffdeltas, field_diffpadvals)
    end
    # global _del = ∇

    # inject sources
    J = inject_sources(source_instances, t)

    # first update E
    # dEdt = (∇ × H - values(E) ⊙ σ - J) ⊘ ϵ

    # tensor subpixel smoothing with staggered gridgp
    # @show typeof.((∇ × H, E ⊙ σ, J))
    global dDdt = (∇ × H - E ⊙ σ - J)
    dEdt = invϵ * dDdt
    E += dEdt * dt

    # then update H
    dBdt = -(∇ × E + values(H) ⊙ m)
    if μ == 1
        dHdt = dBdt
    else
        dHdt = dBdt ⊘ μ
    end
    H += dHdt * dt

    # u = merge(E, H)
    if N == 2
        (Ex=E.Ex, Ey=E.Ey, Hz=H.Hz)
    elseif N == 3
        (Ex=E.Ex, Ey=E.Ey, Ez=E.Ez, Hx=H.Hx, Hy=H.Hy, Hz=H.Hz)
    end
end

# E = apply_field_boundvals(field_boundvals, E; nonzero_only=true)
# """
#     function update(u, p, t, field_boundvals, source_instances; autodiff=true)

# Updates fields. mutating if autodiff is false
# """
# function update(u, p, t, field_deltas, dt, field_boundvals, source_instances; past=false, alg=nothing)
#     # t, field_deltas, dt, field_boundvals, source_instances, autodiff, save_memory = ignore_derivatives() do
#     #     t, field_deltas, dt, field_boundvals, source_instances, autodiff, save_memory
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
#     field_boundvals = ignore_derivatives() do
#         onedge = namedtuple([k => [sum(left, v) sum(right, v)] for (k, v) = pairs(field_boundvals)]) |> cpu
#         1 - onedge
#     end
#     Epads = field_boundvals(r"E.*")
#     ∇ = StaggeredDel(fill(field_deltas, d), field_boundvals, alg)

#     # first update E
#     # dEdt = (∇ × H - E * σ - J) / ϵ
#     # dEdt = invϵ * (∇ × H - E * σ - J)
#     dDdt = (∇ × H - E * σ - J)
#     # dEdt = map(Epads, eachrow(invϵ)) do Epadsi, row
#     dEdt = map(1:size(invϵ, 1)) do i
#         Epadsi = Epads[i]
#         row = invϵ[i, :]
#         # sum(map(Epads, row, dDdt) do Epadsj, a, d
#         sum(map(eachindex(row)) do j
#             Epadsj = Epads[j]
#             a = row[j]
#             d = dDdt[j]
#             l, r = eachcol((Epadsi - Epadsj) / 2)
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

#     # E = apply_field_boundvals(field_boundvals, E; nonzero_only=true)

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

#     # H = apply(field_boundvals, H)
#     # H = apply_field_boundvals(field_boundvals, H; nonzero_only=true)
#     # global asdffsd1=H

#     u = merge(E, H, J0)
#     # _u = merge(E1, H1)
#     # _dudt = merge(dEdt, dHdt)
#     # (_u, _dudt, u)
# end
