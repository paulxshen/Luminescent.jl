
function solve(prob, geometry; autodiff=false, compression=false, ls=nothing, verbose=false, kwargs...)
    @unpack dx, dt, u0, field_padding, geometry_padding, subpixel_averaging, source_instances, monitor_instances, transient_duration, F, polarization, steady_state_duration, d, n = prob

    p = apply_geometry_padding(geometry_padding, geometry)
    # global asdffsd = p
    p = apply_subpixel_averaging(subpixel_averaging, p)
    # return p.ϵxx |> sum

    ignore() do
        # verbose && @show dx, dt, transient_duration, steady_state_duration
        # sz = p |> values |> first |> values |> first |> size
        # points = prod(sz)
        # steps = (transient_duration + steady_state_duration) / dt |> round
        #@debug "" F
        #@debug "size (includes PML): $sz"
        #@debug "$(digitsep(points)) points x $(digitsep(steps)) steps = $(digitsep(points*steps)) point-steps"
    end

    Δ = [transient_duration, steady_state_duration]
    T = cumsum(Δ)
    # N = T[2] / dt |> round
    # milestones = 0:0.1:1.01 |> collect
    # #@debug "simulation started"
    # clock = ignore() do
    #     time()
    # end

    # i = 0
    # _f = ((u,), p, t) -> begin
    #     (update(u, p, t, dx, dt, field_padding, source_instances; autodiff, compression),)
    # end
    # f = (u, t) -> _f(u, p, t)
    if compression
        u, = adjoint_reduce(_f, 0:dt:T[1], (deepcopy(values(u0)),), p, ls)
    else
        u = reduce(0:dt:T[1], init=u0) do u, t
            update(u, p, t, dx, dt, field_padding, source_instances; autodiff, compression)
        end
    end
    return sum.(u) |> sum
    # mf0 = [[
    #     begin
    #         # E = zeros(complex(F), size(m))
    #         # H = zeros(complex(F), size(m))
    #         E = [u[k, m] |> F for k = Enames]
    #         H = [u[k, m] |> F for k = Hnames]
    #         [E, H] * zero(complex(F))
    #     end for λ = m.wavelength_modes |> keys]
    #        for m = monitor_instances
    # ]
    # _f = ((u, mf), p, t) ->
    #     begin
    #         dmf = [[
    #             begin
    #                 E = group(u, :E)
    #                 E = field.((E,), keys(E), (m,))
    #                 H = group(u, :H)
    #                 H = field.((H,), keys(H), (m,))
    #                 [E, H] * cispi(-2t / λ)
    #             end for λ = wavelengths(m)
    #         ] for m = monitor_instances]
    #         mf += dt / Δ[2] * dmf
    #         # if autodiff
    #         # else
    #         #     for (a, b) = zip(mf, dmf)
    #         #         for (a, b) = zip(a, b)
    #         #             for (a, b) = zip(a, b)
    #         #             for (a, b) = zip(a, b)
    #         #                 a .+= dt / Δ[2] * b
    #         #             end
    #         #             end
    #         #         end
    #         #     end
    #         # end
    #         (update(u, p, t, dx, dt, field_padding, source_instances; autodiff, compression), mf,)
    #     end
    # f = (u, t) -> _f(u, p, t)
    ts = T[1]+dt:dt:T[2]
    if compression
        u, mf, = adjoint_reduce(_f, ts, (u, mf0), p, ls)
    else
        u, mf = reduce(ts, init=(u, 0)) do (u, mf), t
            u = update(u, p, t, dx, dt, field_padding, source_instances; autodiff, compression)
            mf += dt / Δ[2] * [[
                begin
                    E = group(u, :E)
                    E = field.((E,), keys(E), (m,))
                    H = group(u, :H)
                    H = field.((H,), keys(H), (m,))
                    [E, H] * cispi(-2t / λ)
                end for λ = wavelengths(m)
            ] for m = monitor_instances]
            (u, mf)
        end
    end

    ls = ignore_derivatives() do
        vcat(extrema.(leaves(u)), [
            begin
                c = maximum(abs.(a))
                (-c, c)
            end for a = leaves(mf)
        ])
    end

    v = map(mf, monitor_instances) do mf, m
        map(mf, values(m.wavelength_modes)) do u, wm
            E, H = u
            if d == 2
                if polarization == :TE
                    Ex, Hy, Ez = invreframe(frame(m), vcat(E, H))
                    Ex += 0sum(Ez[1:2])
                    Hy += 0sum(Ez[1:2])
                    mode = (; Ex, Hy, Ez)
                else
                    Hx, Ey, Hz = invreframe(frame(m), vcat(H, E))
                    Hx += 0real(Hz[1])
                    Ey += 0real(Hz[1])
                    mode = (; Hx, Ey, Hz)
                end
            elseif d == 3
                Ex, Ey, Ez, = invreframe(frame(m), E)
                Hx, Hy, Hz = invreframe(frame(m), H)
                mode = (; Ex, Ey, Ez, Hx, Hy, Hz)
            end
            # mode = keepxy(mode)
            c = mode_decomp.(wm, (mode,), dx)
            fp = [abs(v[1])^2 for v = c]
            rp = [abs(v[2])^2 for v = c]
            mode, c, fp, rp
        end
    end
    modes = [[v[1] for v = v] for v in v]
    mode_coeffs = [[v[2] for v = v] for v in v]
    forward_mode_powers = [[v[3] for v = v] for v in v]
    reverse_mode_powers = [[v[4] for v = v] for v in v]

    return (; fields=u, geometry, modes, mode_coeffs, forward_mode_powers, reverse_mode_powers, ls)
end