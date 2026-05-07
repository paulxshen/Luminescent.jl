bboxed(a, Ibbox) =
    fmap(a) do a
        a(Ibbox...)
    end

Zygote._unapply(xs::Tuple, t::Tuple{}) = xs, ()
Zygote._unapply(::Tuple{}, ::Tuple{}) = (), ()


function stepper(u, p, udft; i, t, dt, umax=nothing, Ibbox, source_instances, monitor_instances, saveat, checkat, ϵ, weights, verbosity, FIELDS, FIELDS_SPECIAL, ipeak, monitorat, ∇ₕ, ∇ₑ, ks, fs, _dt, ts, energy_decay_threshold)
    @nograd i, t, dt, source_instances, monitor_instances, saveat, ϵ, weights, verbosity, FIELDS, ipeak, monitorat, ∇ₕ, ∇ₑ, ks, fs, _dt, ts, energy_decay_threshold

    ignore_derivatives() do
        checked = false
        if i % checkat == 0
            checked = true
            # print("")
            if any(isnan, u(1)(1))
                error("instability. try lower `relative_courant`")
                # exit(3)
            end
            if !isnothing(energy_decay_threshold)
                E_ = energy(u, ϵ, weights, Ibbox)
                if E_ < energy_decay_threshold * Emax
                    global stop = true
                end
            end
            check_ping()
        end
        if i == ipeak
            checked = true
            global Emax = 2energy(u, ϵ, weights, Ibbox)
        end
        if !isnothing(umax) && checked
            for k1 = (:E, :H)
                for (k, v) = pairs(u(k1))
                    umax[k] = max(get(umax, k, 0), maximum(abs, v))
                end
            end
        end
    end

    if (!isnothing(saveat) && i % saveat == 1) || i == length(ts) || i == ipeak || stop
        ignore_derivatives() do
            name = "$(round3(t)).npz"
            FS = []
            if i == ipeak
                push!(FS, FIELDS_SPECIAL)
                # name = "peak$SEP$name"
            elseif i == length(ts) || stop
                push!(FS, FIELDS_SPECIAL)
                # name = "final$SEP$name"
            end
            if !isnothing(saveat) && i % saveat == 1
                push!(FS, FIELDS)
            end

            if verbosity > 1
                E_ = energy(u, ϵ, weights, Ibbox)
                tp = timepassed()
                if i > ipeak
                    s = "energy $(E_/Emax|>round4), "
                else
                    s = ""
                end
                println("period $(t|>disp), $s$(tp/(t-_t)|>disp) seconds/period")
                global _t = t
            end
            _u = merge(u[:E], u[:H])
            for x = FS
                try
                    npzwrite(joinpath(x, name), Dict(string.(keys(_u)) .=> cpu(_values(_u)) |> FF))
                catch e
                    println(e)
                end
            end
        end
    end
    if stop
        println("stop on field decay")
        return
    end

    if i % monitorat == 0
        um = map(monitor_instances) do m
            @nograd m
            namedtuple(ks .=> field.((u,), ks, (m,)))
        end
        c = _dt * cispi.(-2t * fs)
        if AD()
            udft += map(um) do v
                _c = reshape(c, ones(Int, ndims(v(1)))..., length(c))
                vmap(v) do v
                    _c .* v
                end
            end
        else
            map(um, udft) do v, udft
                _c = reshape(c, ones(Int, ndims(v(1)))..., length(c))
                map(ks) do k
                    # global _a = udft, k, _c, v
                    udft[k] .= udft[k] + _c .* v[k]
                end
            end
        end
        # fluxdft += map(um) do u
        #     S = stack(poynting(u))
        #     _c = reshape(c, ones(Int, ndims(S))..., length(c))
        #     _c * S
        # end
    end

    u = update(u, p, t; dt, ∇ₕ, ∇ₑ, source_instances)
    u, p, udft
end

function foo(models; prob)
    ENV["timestamp"] = prob[:timestamp]
    println("""

    $BREAK
    running simulation...
    """)


    @unpack approx_2D_mode, dt, u0, geometry, source_instances, path, monitor_instances, centgeom, design_instances, tmax, meshprops, energy_decay_threshold, grid, saveat, checkat, verbosity, relative_courant, λs, fs, hasPEC, gradckptat, array, umax, time_extrapolation = prob

    # ENV["timestamp"] = prob[:timestamp]

    @nograd dt, u0, source_instances, monitor_instances, tmax, grid
    @unpack ∇ₕ, ∇ₑ, offsets, F, N, weights, Ibbox, Ibbox1 = grid

    dt = 1 / ceil(1 / dt) |> FF
    n = hasPEC ? 64 : 32
    monitorat = @ignore_derivatives max(1, floor(Int, 1 / n / dt))

    if !isnothing(saveat)
        saveat = max(1, round(Int, saveat / dt))
    end

    v = call.(models)
    v = [d.partial_etch ? min.(v[i], v[i-1]) : v[i] for (i, d) in enumerate(design_instances)]
    for (c, m) = zip(design_instances, v)
        println("applying design: $(c.name)")
        geometry = apply_design!(c, geometry, centgeom, offsets, m, meshprops)
    end
    p = geometry

    FIELDS = joinpath(path, "fields")
    FIELDS_SPECIAL = joinpath(path, "fields_special")
    rm(FIELDS, force=true, recursive=true)
    rm(FIELDS_SPECIAL, force=true, recursive=true)
    mkpath(FIELDS)
    mkpath(FIELDS_SPECIAL)

    ignore_derivatives() do
        d = cpu(centgeom)
        npzwrite(joinpath(path, "g.npz"), Dict(string.(keys(d)) .=> _values(d) |> FF))
        #  CairoMakie.save(joinpath(path, "geometry.png"), 
    end

    timepassed()
    ts = range(0, tmax - dt, round(Int, tmax / dt))
    ipeak = searchsortedfirst(ts, maximum(duration.(source_instances)) / 2)
    @nograd ts

    verbosity > 1 && println("accumulating dft fields...")
    u = u0
    ks = @ignore_derivatives [keys(u.E)..., keys(u.H)...]

    udft = @ignore_derivatives map(monitor_instances) do m
        if isa(m, PlaneMonitorInstance)
            NamedTuple([k => zeros(complex(FF), length.(m.Is[k])..., length(fs)) |> array for k = ks])
        else
            NamedTuple([k => zeros(complex(FF), length(m.Is[k][1]), length(fs)) |> array for k = ks])
        end
    end
    # udft = 0
    ϵ = geometry(:ϵ)[1]
    _dt = monitorat * dt
    global _t = 0
    global stop = false
    global Emax = 0f0

    # @debug gradckptat#checkp
    kw = (; ts, source_instances, monitor_instances, saveat, checkat, ϵ, weights, verbosity, FIELDS, FIELDS_SPECIAL, ipeak, monitorat, ∇ₕ, ∇ₑ, ks, fs, dt, _dt, energy_decay_threshold, Ibbox=Ibbox1)
    tf = tmax

    if !AD()
        ignore_derivatives() do
            #     cache = GPUArrays.AllocCache()
            for (i, t) = enumerate(ts)
                # GPUArrays.@cached cache begin
                _u = stepper(u, p, udft; i, t, umax, kw...)
                stop && break
                u, p, udft = _u
                tf = t
                # end
            end
        end
    else
        n = length(ts) ÷ gradckptat
        U = A = 0
        ignore_derivatives() do
            U = [fmap(similar, u) for _ in 1:n]
            # U = nothing
            A = [fmap(similar, udft) for _ in 1:n]
            l = prod(size(u(:Ex)))
            @debug "malloc" 24f-9l * n

        end
        u, p, udft = adjoint_reduce(stepper, u, p, udft; U, A, gradckptat, umax, array, kw...)
    end

    if time_extrapolation
        deconv = @ignore_derivatives Deconv(fs |> Array |> sort, tf)
    else
        deconv = identity
    end

    c = cispi.(-dt * reverse(fs))
    udft = map(udft) do d
        kvmap(d) do k, v
            _c = reshape(c, fill(1, ndims(v) - 1)..., length(c))
            k => startswith(string(k), "H") ? _c .* v : v
        end
    end
    udft = map(udft, monitor_instances) do d, m
        @nograd m
        d = localframe(d, m; approx_2D_mode)
    end
    fields = map(udft) do d
        n = size(d(1))[end]
        map(1:n) do i
            vmap(d) do a
                selectdim(a, ndims(a), i)
            end
        end
    end

    # flux = map( monitor_instances) do a, m
    #     if isa(m, PlaneMonitorInstance)
    #         @show size(a)
    #         a = map(m.local_plane_Is) do I
    #             getindexf(a, I..., :, :)
    #         end
    #         @show size(a)
    #         local N = ndims(a)
    #         nf = size(a)[end]
    #         sz = size(a)[1:end-2]
    #         a = permutedims(a, (N - 1, (1:N-2)..., N))
    #         a = reshape(a, size(a, 1), :)
    #         a = eachcol(m.frame)[end]' * a
    #         a = reshape(a, sz..., nf)
    #         sum(m.dA .* a, dims=1:ndims(a)-1)
    #     else
    #     end
    # end
    _flux = map(udft, monitor_instances) do u, m
        a = m.dA .* real(conj(u.Ex) .* u.Hy - conj(u(:Ey, 0)) .* u(:Hx, 0))
        vec(sum(a, dims=1:ndims(a)-1)) |> Array
    end
    _waves = map(fields, monitor_instances) do us, m
        @nograd m
        @unpack dA = m
        isnothing(m.λmodes) && return []
        map(eachcol(m.λmodes), eachcol(m._λmodes)) do modes, _modes
            [inner.(modes, us, (dA,)), inner.(_modes, us, (dA,))]
        end
    end
    _waves = fmap(deconv, _waves)

    profiles = map(fields, monitor_instances) do ds, m
        map(ds) do d
            _inner.((abs2,), m.modes, (d,)), _inner.((abs2,), m._modes, (d,))
        end #|> stack
    end

    flux = OrderedDict()
    flux_net = OrderedDict()
    flux_decomp = OrderedDict()
    flux_radiation = OrderedDict()
    waves = OrderedDict()
    radiation = OrderedDict()
    radiation_isotropic = OrderedDict()
    pm = @ignore_derivatives "+", "-"
    HLVR = @ignore_derivatives "H", "L", "V", "R"
    # sensors = OrderedDict(["V" => OrderedDict(), "I" => OrderedDict(), "Z" => OrderedDict()])
    V = OrderedDict()
    I = OrderedDict()
    Z = OrderedDict()
    for (mon, w, f, p, λu,) = zip(monitor_instances, _waves, _flux, profiles, fields)
        @unpack name, A, dA = mon
        if mon isa PlaneMonitorInstance
            @unpack voltage_sample_Is, voltage_dl, current_sample_Iss, current_dls, λZ = mon

            if !isnothing(voltage_sample_Is)
                V[name] = map(λu) do d
                    do_line_integral([d(:Ex), d(:Ey)], voltage_sample_Is, voltage_dl)
                end
            end
            if !isnothing(current_sample_Iss)
                I[name] = map(λu) do d
                    sum(zip(current_sample_Iss, current_dls)) do (pts, dl)
                        do_line_integral([d(:Hx), d(:Hy)], pts, dl)
                    end
                end
            end
            if !isnothing(λZ)
                # Z = V ./ I
                # sensors["Z"][name] = Z
                Z[name] = λZ
                v = map(V[name], I[name], λZ) do v, i, z
                    [1 1; 1/z -1/z] \ [v, i]
                end
                waves["$name+"] = deconv(first.(v))
                waves["$name-"] = deconv(last.(v))
            end
        end

        name = string(name)
        flux[name] = f
        for (i, dir) = enumerate(pm)
            v = getindex.(p, i)
            for (i, pol) = enumerate(HLVR)
                radiation["$name@$pol$dir"] = getindex.(v, i)
            end
            flux_radiation["$name$dir"] = map(v) do v
                sum(dot.(v, (dA,)))
            end
        end

        flux_radiation[name] = flux_radiation["$name+"] - flux_radiation["$name-"]

        for dir = pm
            radiation["$name$dir"] = sum(HLVR) do pol
                radiation["$name@$pol$dir"]
            end
        end
        radiation[name] = radiation["$name+"] - radiation["$name-"]

        for pol = HLVR
            radiation_isotropic["$name@$pol"] = (radiation["$name@$pol+"] - radiation["$name@$pol-"]) ./ flux_radiation[name] * A
        end
        radiation_isotropic[name] = sum(HLVR) do pol
            radiation_isotropic["$name@$pol"]
        end

        for (i, dir) = enumerate(pm)
            if isa(mon, PlaneMonitorInstance)
                flux_decomp["$name$dir"] = 0
                for (m, v) = enumerate(w)
                    waves["$(name)@$(m-1)$dir"] = v[i]
                    flux_decomp["$name$dir"] += abs2.(v[i])
                end
                flux["$name$dir"] = flux_decomp["$name$dir"]
            elseif isa(mon, SphereMonitorInstance)
                flux["$name$dir"] = flux_radiation["$name$dir"]
            end
        end
    end

    # AD() && return 
    # (; waves, flux, umax)
    #nc
    radiation = vmap(stack, radiation)
    radiation_isotropic = vmap(stack, radiation_isotropic)

    println("\nsimulation completed.")
    (; waves, fields, flux, profiles, flux_decomp, flux_radiation, radiation, radiation_isotropic, umax, V, I, Z)
end


struct Solution
    u
    p
    udft
    P
    ap
    am
end
@functor Solution

function (s::Solution)(k, m, w=1, mn=0)
    @unpack u, P, udft, ap, am = s
    if k == "a+"
        return s.ap[m][w][mn+1]
    elseif k == "a-"
        return s.am[m][w][mn+1]
    elseif k == "P_TE"
        return flux(udft[m][w], :TE)
    elseif k == "P_TM"
        return flux(udft[m][w], :TM)
    elseif k == "P"
        return P[m][w]
    elseif k == "udft"
        return udft[m][w]
    end
end

