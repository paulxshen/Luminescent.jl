
function solve(prob, ; autodiff=true, history=nothing, save=false)
    @unpack dx, dt, u0, field_padding, source_instances, monitor_instances, transient_duration, steady_state_duration, geometry, wavelengths, modes = prob

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

    u, fp, = reduce(T[1]+dt:dt:T[2], init=(u, 0, 0)) do (u, fp,), t
        save && push!(h, u)
        fp_ = [
            begin
                E = [u[k, m] for k = Enames]
                H = [u[k, m] for k = Hnames]
                [E, H]
            end
            for m = monitor_instances
        ]
        # tp_ = map(fp_, monitor_instances) do (E, H), m
        #     mean(sum((E × H) .* normal(m))) * area(m)
        # end

        (
            _update(u, p, t, dx, dt, field_padding, source_instances),
            [fp + dt / Δ[2] * fp_ * cispi(-2t / λ) for λ = wavelengths],
            # tp + dt / Δ[2] * tp_,
        )
    end

    v = map(fp, modes) do fp, mode
        map(fp, monitor_instances) do (E, H), m
            if polarization == :TE
                # global Ex, Hy
                Ex, Hy, Ez = invreframe(frame(m), vcat(E, H))
                Ex += 0real(Ez[1])
                Hy += 0real(Ez[1])
            else
            end
            mode_decomp(mode, (; Ex, Hy))
            # @show abs(ap)^2, abs(am)^2
            # real((Ex ⋅ nfp.Ex) ⋅ (Hy ⋅ nfp.Hy)) * nmp + 0real(Ez[1])
            # abs(Ex ⋅ nfp.Ex) ⋅ abs(Hy ⋅ nfp.Hy) * nmp + 0real(Ez[1])
        end
    end
    stack(stack.(v))
    # getindex.(values(v), 1), getindex.(values(v), 2)
end