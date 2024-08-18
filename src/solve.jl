
function solve(prob; autodiff=false, compression=false, lims=nothing, history=nothing, comprehensive=true, verbose=false, plotpath=nothing, _time=true,
    cpu=identity, gpu=identity, kwargs...)
    @unpack dx, dt, u0, field_padding, geometry_padding, subpixel_averaging, source_instances, geometry, monitor_instances, transient_duration, F, polarization, steady_state_duration, d, n = prob


    p = geometry
    _gpu = isa(first(values(p)), AbstractGPUArray)
    p = apply(geometry_padding, p)
    p = apply(subpixel_averaging, p)

    ignore() do
        # verbose && @show dx, dt, transient_duration, steady_state_duration
        sz = p |> values |> first |> values |> first |> size
        points = prod(sz)
        steps = (transient_duration + steady_state_duration) / dt |> round
        #@debug "" F
        #@debug "size (includes PML): $sz"
        #@debug "$(digitsep(points)) points x $(digitsep(steps)) steps = $(digitsep(points*steps)) point-steps"
    end

    Δ = [transient_duration, steady_state_duration]
    T = cumsum(Δ)
    N = T[2] / dt |> round
    Enames = keys(u0.E)
    Hnames = keys(u0.H)
    # extrema(abs.(prob.source_instances[1]._g.Jy))

    # if save
    milestones = 0:0.1:1.01 |> collect
    #@debug "simulation started"
    clock = ignore() do
        time()
    end

    _reduce = (autodiff && compression) ? adjoint_reduce : reduce
    i = 0
    u = _reduce(0:dt:T[1], (deepcopy(u0),), lims) do (u,), t
        verbose && ignore() do
            hasnan(u) && error("nan detected. instability. aborting")
            if !isempty(milestones) && (i + 1) / N > milestones[1]
                println("$(milestones[1]*100)% done $(time()-clock)s since start")
                deleteat!(milestones, 1)
            end
            i += 1
        end
        (update(u, p, t, dx, dt, field_padding, source_instances; autodiff),)
    end

    # u, mode_fields, total_powers = reduce(T[1]+dt:dt:T[2], init=(u, 0, 0)) do (u, mode_fields, tp), t
    u, mode_fields, = _reduce(T[1]+dt:dt:T[2], (u, 0), lims) do (u, mode_fields), t
        verbose && ignore() do
            if !isempty(milestones) && (i + 1) / N > milestones[1]
                hasnan(u) && error("nan detected. instability. aborting")
                println("$(milestones[1]*100)% done $(time()-clock)s since start")
                deleteat!(milestones, 1)
            end
            i += 1
        end

        mode_fields_ = [
            [
                begin
                    E = [u[k, m] |> F for k = Enames]
                    H = [u[k, m] |> F for k = Hnames]
                    [E, H] * cispi(-2t / λ)
                end for λ = m.wavelength_modes |> keys]
            for m = monitor_instances
        ]
        if autodiff
            mode_fields += dt / Δ[2] * mode_fields_
        else
            for (a, b) = zip(mode_fields, mode_fields_)
                for (a, b) = zip(a, b)
                    for (a, b) = zip(a, b)
                        a .+= dt / Δ[2] * b
                    end
                end
            end
        end
        # tp_ = map(monitor_instances) do m
        #     E = [u[k, m] for k = Enames]
        #     H = [u[k, m] for k = Hnames]
        #     mean(sum((E × H) .* normal(m))) * area(m)
        # end

        (
            update(u, p, t, dx, dt, field_padding, source_instances; autodiff),
            mode_fields,
            # tp + dt / Δ[2] * tp_,
        )
    end
    #@debug "simulation ended"
    # ignore() do
    #     println("simulation took: ", time() - clock, "s")
    # end

    lss = ignore_derivatives() do
        extrema.(leaves(u)), [
            begin
                c = maximum(abs.(a))
                (-c, c)
            end for a = leaves(mode_fields)
        ]
    end

    v = map(mode_fields, monitor_instances) do mf, m
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

    fields = u
    # @show forward_mode_powers, reverse_mode_powers, total_powers
    # return (; fields, geometry, modes, mode_coeffs, forward_mode_powers, reverse_mode_powers, total_powers)
    return (; fields, geometry, modes, mode_coeffs, forward_mode_powers, reverse_mode_powers, lss)
end