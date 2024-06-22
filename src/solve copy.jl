
function solve(prob, model=nothing; autodiff=true, history=nothing, comprehensive=true, verbose=false, plotpath=nothing, kwargs...)
    @unpack dx, dt, u0, field_padding, geometry_padding, subpixel_averaging, source_instances, monitor_instances, transient_duration, steady_state_duration, d = prob
    # if isnothing(geometry)
    # end

    p = prob.geometry
    mask = model[1]()
    f = design_config.fill.ϵ
    v = design_config.void.ϵ
    # for (f, v) = zip(design_config.fill, design_config.void)
    a = (v .* (1 - mask) + f .* mask)
    ϵ = place(p[:ϵ], ((designs[1].bbox[1] - origin) / dx) + 1, a; replace=true) |> F

    # p = cpu(p)
    p = apply(geometry_padding, (; ϵ, μ=p.μ, σ=p.σ, σm=p.σm))
    p = apply(subpixel_averaging, p)
    geometry = p
    # p = gpu(p)

    Δ = [transient_duration, steady_state_duration]
    T = cumsum(Δ)
    Enames = keys(u0.E)
    Hnames = keys(u0.H)
    # extrema(abs.(prob.source_instances[1]._g.Jy))

    # if save
    _update = !autodiff && !save ? maxwell_update! : maxwell_update
    h = []

    # run simulation
    u = reduce(0:dt:T[1], init=deepcopy(u0)) do u, t
        # save && push!(h, u)
        _update(u, p, t, dx, dt, field_padding, source_instances;)
    end

    u, mode_fields, total_powers = reduce(T[1]+dt:dt:T[2], init=(u, 0, 0)) do (u, mode_fields, tp), t
        # save && push!(h, u)
        mode_fields_ = [
            [begin
                E = dict([k => u[k, m] for k = Enames])
                H = dict([k => u[k, m] for k = Hnames])
                (; E, H)
            end * cispi(-2t / λ) for λ = m.wavelength_modes |> keys]
            for m = monitor_instances
        ]
        tp_ = map(monitor_instances) do m
            E = [u[k, m] for k = Enames]
            H = [u[k, m] for k = Hnames]
            mean(sum((E × H) .* normal(m))) * area(m)
        end

        (
            _update(u, p, t, dx, dt, field_padding, source_instances),
            mode_fields + dt / Δ[2] * mode_fields_,
            tp + dt / Δ[2] * tp_,
        )
    end

    v = map(mode_fields, monitor_instances) do v, m
        map(v, m.wavelength_modes |> keys) do u, λ
            E = values(u.E)
            H = values(u.H)
            if d == 2
                if polarization == :TE
                    # global E, H
                    Ex, Hy, Ez = invreframe(frame(m), vcat(E, H))
                    Ex += 0real(Ez[1])
                    Hy += 0real(Ez[1])
                    _mode = (; Ex, Hy, Ez)
                else
                    Hx, Ey, Hz = invreframe(frame(m), vcat(H, E))
                    Hx += 0real(Hz[1])
                    Ey += 0real(Hz[1])
                    _mode = (; Hx, Ey, Hz)
                end
            elseif d == 3
                Ex, Ey, Ez, = invreframe(frame(m), E)
                Hx, Hy, Hz = invreframe(frame(m), H)
                _mode = (; Ex, Ey, Ez, Hx, Hy, Hz)
            end
            _mode = keepxy(_mode)
            c = mode_decomp.(m.wavelength_modes[λ], (_mode,), dx)
            fp = [abs(v[1])^2 for v = c]
            rp = [abs(v[2])^2 for v = c]
            _mode, c, fp, rp
        end
    end
    modes = [[v[1] for v = v] for v in v]
    mode_coeffs = [[v[2] for v = v] for v in v]
    forward_mode_powers = [[v[3] for v = v] for v in v]
    reverse_mode_powers = [[v[4] for v = v] for v in v]

    fields = u
    if verbose
        @show forward_mode_powers, reverse_mode_powers, total_powers
    end
    if !isnothing(plotpath)

    end
    if comprehensive
        return (; fields, geometry, modes, mode_coeffs, forward_mode_powers, reverse_mode_powers, total_powers)
    end
    forward_mode_powers
end