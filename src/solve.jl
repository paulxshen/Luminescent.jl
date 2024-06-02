
function solve(prob, ; autodiff=true, history=nothing, save=false, comprehensive=false, verbose=false)
    @unpack dx, dt, u0, field_padding, source_instances, monitor_instances, transient_duration, steady_state_duration, geometry = prob

    p = geometry
    p = apply(geometry_padding; p)
    p = apply(subpixel_averaging, p)

    Δ = [transient_duration, steady_state_duration]
    T = cumsum(Δ)
    Enames = keys(u0.E)
    Hnames = keys(u0.H)

    # if save
    _update = !autodiff && !save ? maxwell_update! : maxwell_update
    h = []

    # run simulation
    u = reduce(0:dt:T[1], init=deepcopy(u0)) do u, t
        save && push!(h, u)
        _update(u, p, t, dx, dt, field_padding, source_instances;)
    end

    u, fields, total_powers = reduce(T[1]+dt:dt:T[2], init=(u, 0, 0)) do (u, fields, tp), t
        save && push!(h, u)
        fields_ = [
            [begin
                E = [u[k, m] for k = Enames]
                H = [u[k, m] for k = Hnames]
                [E, H]
            end
             *
             cispi(-2t / λ) for λ = m.wavelengths]
            for m = monitor_instances
        ]
        tp_ = map(fields_, monitor_instances) do (E, H), m
            mean(sum((E × H) .* normal(m))) * area(m)
        end

        (
            _update(u, p, t, dx, dt, field_padding, source_instances),
            fields + dt / Δ[2] * fields_,
            tp + dt / Δ[2] * tp_,
        )
    end

    mode_powers = map(fields, monitor_instances) do v, m
        map(v, m.wavelengths, m.modes) do (E, H), λ, mode
            if d == 2
                if polarization == :TE
                    # global Ex, Hy
                    Ex, Hy, Ez = invreframe(frame(m), vcat(E, H))
                    Ex += 0real(Ez[1])
                    Hy += 0real(Ez[1])
                    mode_decomp(mode, (; Ex, Hy))
                else
                    Hx, Ey, Hz = invreframe(frame(m), vcat(H, E))
                    Hx += 0real(Hz[1])
                    Ey += 0real(Hz[1])
                    mode_decomp(mode, (; Hx, Ey))
                end
            else
                Ex, Ey, Ez, = invreframe(frame(m), E)
                Hx, Hy, Hz = invreframe(frame(m), H)
                mode_decomp(mode, (; Ex, Hy, Hx, Ey))
            end
            # @show abs(ap)^2, abs(am)^2
            # real((Ex ⋅ nfields.Ex) ⋅ (Hy ⋅ nfields.Hy)) * nmp + 0real(Ez[1])
            # abs(Ex ⋅ nfields.Ex) ⋅ abs(Hy ⋅ nfields.Hy) * nmp + 0real(Ez[1])
        end
    end
    forward_mode_powers = [[v[1] for v = v] for v in mode_powers]
    reverse_mode_powers = [[v[2] for v = v] for v in mode_powers]

    if verbose
        @show forward_mode_powers, reverse_mode_powers, total_powers
    end
    if comprehensive
        return (; fields, forward_mode_powers, reverse_mode_powers, total_powers)
    end
    forward_mode_powers
end