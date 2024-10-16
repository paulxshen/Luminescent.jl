function f1(((u,), p, (dx, dt, field_padvals, source_instances)), t)
    u = update(u, p, t, dx, dt, field_padvals, source_instances)
    ((u,), p, (dx, dt, field_padvals, source_instances))
end

function f2(((u, mf), p, (dx, dt, field_padvals, source_instances), (T, monitor_instances)), t)
    u = update(u, p, t, dx, dt, field_padvals, source_instances;)
    mf += dt / T * [[
        begin
            E = u(r"E.*")
            E = field.((E,), keys(E), (m,))

            H = u(r"H.*")
            H = field.((H,), keys(H), (m,))
            [E, H] * cispi(-2t / λ)
        end for λ = wavelengths(m)
    ] for m = monitor_instances]
    ((u, mf), p, (dx, dt, field_padvals, source_instances), (T, monitor_instances))
end

function solve(prob, ;
    save_memory=false, ulims=(-3, 3), framerate=0, path="",
    kwargs...)
    @unpack dx, dt, u0, geometry, _geometry, field_padvals, geometry_padvals, geometry_padamts, fieldlims, source_instances, monitor_instances, transient_duration, F, polarization, steady_state_duration, N, sz, ratio = prob
    ratio = 1
    # transient_duration = steady_state_duration = 0.1
    u0, dx, dt, field_padvals, geometry_padvals, geometry_padamts, fieldlims, source_instances, monitor_instances, transient_duration, F, polarization, steady_state_duration, N, sz, ratio = ignore_derivatives() do
        u0, dx, dt, field_padvals, geometry_padvals, geometry_padamts, fieldlims, source_instances, monitor_instances, transient_duration, F, polarization, steady_state_duration, N, sz, ratio
    end
    fieldlims = cpu(fieldlims)

    p = pad_geometry(geometry, geometry_padvals, geometry_padamts, ratio)
    # global a1 = p
    # _p = pad_geometry(_geometry, geometry_padvals, ratio)

    p = apply_subpixel_averaging(p, fieldlims)
    # _p = apply_subpixel_averaging(_p, fieldlims)
    global a2 = p

    # global _ϵ = ϵ
    # invϵ = tensorinv(_p.ϵ, ratio, fieldlims)

    # p = merge(p, (; invϵ))
    # p = merge(p, _p)
    Δ = [transient_duration, steady_state_duration]
    T = cumsum(Δ)
    us0 = (u0,)
    init = (us0, p, (dx, dt, field_padvals, source_instances))

    if save_memory
        (u,), = adjoint_reduce(f1, 0:dt:T[1], init, ulims)
    else
        (u,), = reduce(0:dt:T[1]; init) do us, t
            ignore() do
                if framerate > 0 && t > 0
                    if t % (1 / framerate) < dt
                        (u,), p, = us
                        a = u.Hz
                        g = p.ϵxx

                        _path = joinpath(path, "temp")
                        mkpath(_path)
                        CairoMakie.save(joinpath(_path, "$t.png"), quickie(a, g; monitor_instances, source_instances, ulims),)
                        # quickie(a, g; monitor_instances, source_instances)
                    end
                end
            end
            f1(us, t)
        end
    end
    ts = T[1]+dt:dt:T[2]+F(0.001)
    init = ((u, 0), p, (dx, dt, field_padvals, source_instances), (Δ[2], monitor_instances))

    if save_memory
        (u, mf), = adjoint_reduce(f2, ts, init, ulims)
    else
        (u, mf), = reduce(f2, ts; init)
    end

    ulims = ignore_derivatives() do
        map(vcat(extrema.(leaves(u)), [
            begin
                c = maximum(abs.(a))
                (-c, c)
            end for a = leaves(mf)
        ])) do l
            l[1] == l[2] ? convert.(F, (-1, 1)) : l
        end
    end

    v = map(mf, monitor_instances) do mf, m
        map(mf, wavelengths(m)) do u, λ
            E, H = u
            if N == 2
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
            elseif N == 3
                Ex, Ey, Ez, = invreframe(frame(m), E)
                Hx, Hy, Hz = invreframe(frame(m), H)
                mode = (; Ex, Ey, Ez, Hx, Hy, Hz)
            end
            # mode = keepxy(mode)
            fp = rp = c = nothing
            if !isnothing(m.wavelength_modes)
                wm = m.wavelength_modes[λ]
                c = mode_decomp.(wm, (mode,), dx)
                fp = [abs(v[1])^2 for v = c]
                rp = [abs(v[2])^2 for v = c]
            end

            mode, c
        end
    end
    um = [[v[1] for v = v] for v in v]
    ap = [[getindex.(v[2], 1) for v = v] for v in v]
    am = [[getindex.(v[2], 2) for v = v] for v in v]
    return Solution(u, p, ulims, um, ap, am)
end

struct Solution
    u
    p
    ulims
    um
    ap
    am
end
@functor Solution

function (s::Solution)(k, m, w=1, mn=0)
    @unpack u, ulims, um, ap, am = s
    if k == "a+"
        return s.ap[m][w][mn+1]
    elseif k == "a-"
        return s.am[m][w][mn+1]
    elseif k == "P_TE"
        return flux(um[m][w], :TE)
    elseif k == "P_TM"
        return flux(um[m][w], :TM)
    elseif k == "P"
        return flux(um[m][w])
    elseif k == "um"
        return um[m][w]
    end
end
