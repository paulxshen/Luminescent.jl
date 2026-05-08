function picrun(path; array=Array, kwargs...)
    # $(abspath(path)) 
    println("""
    $BREAK
    backend loading simulation folder ...
    """)
    _versioninfo()

    Random.seed!(1)
    PROB = joinpath(path, "problem.json")
    SOL = joinpath(path, "solution.json")

    prob = readjsonnp(PROB)
    for (k, v) = kwargs
        prob[k] = v
    end


    @unpack name, N, approx_2D_mode, dtype, boundaries, wavelength, wavelengths, monitors, sources, relative_courant, relative_pml_depths, z, magic, layer_stack, material_library, tmax, energy_decay_threshold, path_length_multiple, bbox, nres, dx, dy, dz, modes, saveat, gradient_checkpoint, ordering, secret, subpixel_smoothing, subpixel_smoothing_sampling_distance, designs, load_saved_designs, targets, time_extrapolation, timestamp, study_type = prob
    pro = hash(secret) == 0x31634f070f133a63
    # pro = true
    designs = map(enumerate(designs)) do (i, c)
        @unpack name, uniform_along, symmetries, swaps, lmin, symmetries, partial_etch = c
        swaps = kvmap(swaps) do k, v
            Symbol(togreek(k)) => v
        end
        Design(name, swaps, stack(c(:bbox)), lmin; uniform_along, symmetries, partial_etch)
    end



    if study_type == "inverse_design"
        println("inverse design not available on community cloud. Please upgrade to a private cloud.")
        # error()

        verbosity = 1
        @unpack optimizer, pixel_size = prob
        @unpack iters, stoploss, lowloss, momentum, contrast, ckptat, gradckptat = optimizer

        topopt_kwargs = (; targets)
    else
        verbosity = 2
        gradckptat = nothing
        topopt_kwargs = (;)
    end
    global _prob0 = prob
    λ = FF(wavelength)
    λs = FF.(wavelengths)
    bbox = FF.(stack(bbox))

    layer_stack = sort(collect(pairs(layer_stack)), by=kv -> kv[2](:mesh_order)) |> OrderedDict
    ks = keys(layer_stack)
    fns = readdir(joinpath(path, "geometry"), join=true)
    sort!(fns, by=x -> parse(Int, split(basename(x), SEP)[1]))
    @debug (; material_library)
    material_library = vmap(material_library) do d
        kvmap(d) do k, v
            togreek(Symbol(k)), FF(v)
        end
    end
    bg = material_library("bg")
    geometry = map(fns) do f
        m = getfield(GeoIO.load(f; numtype=Float32), :domain)
        v = split(basename(f)[1:end-4], SEP)
        d = material_library(v[2]) |> pairs |> OrderedDict
        m, d
    end

    if N == 2
        bbox = bbox[1:2, :]
    end
    if study_type == "inverse_design"
        ENV["AD"] = 1
        if isfile(SOL)
            sol = open(SOL, "r") do f
                JSON.parse(f)
            end
        else
            sol = nothing
        end
    else
        ENV["AD"] = 0
    end

    sources = map(sources) do source
        @unpack name, origin, frame, start, stop, modenums, wavelength, bandwidth, duration = source
        origin = origin[1:N]
        start = start[1:N-1]
        stop = stop[1:N-1]

        if N == 2
            frame = [frame[1:2, 1] frame[1:2, 3]]
        end

        Source(wavelength, bandwidth, duration, origin, start, stop, frame; modenums, name, label="source_$name")
    end

    monitors = map(monitors) do m
        label = m(:name)
        if string(m(:type)) == "plane"
            @unpack name, origin, frame, start, stop = m
            origin = origin[1:N]
            start = start[1:N-1]
            stop = stop[1:N-1]

            if N == 2
                frame = [frame[1:2, 1] frame[1:2, 3]]
            end

            PlaneMonitor(name, origin, start, stop, frame; label)
        elseif m(:type) == "sphere"
            @unpack name, origin, radius, frame, θmin, θmax, φmin, φmax, angres = m

            SphereMonitor(name, origin, radius, frame; label, θmin, θmax, φmin, φmax, angres)
        else
            error("unknown monitor type $(m(:type))")
        end
    end

    F = FF
    println()
    global prob = setup(λ, λs, bbox, nres, geometry, boundaries, sources, monitors, modes;
        relative_pml_depths, approx_2D_mode, z,
        saveat, verbosity, relative_courant,
        F, tmax, energy_decay_threshold, path_length_multiple, path,
        subpixel_smoothing, name,
        bg, time_extrapolation, timestamp,
        designs, targets, gradckptat, array
    )
    @unpack design_instances = prob
    models = getfield.(design_instances, :model)

    if load_saved_designs
        println("loading saved designs...")
        for (m, d) = zip(models, design_instances)
            m.p .= resize(npzread(joinpath(path, "designs", "params_$(d.name).npy")), size(m.p))
        end
    end


    if study_type == "sparams"
        ENV["AD"] = 0
        global sol = foo(models; prob)
        @unpack monitor_instances = prob
        @unpack flux, fields, waves, flux_decomp, flux_radiation, radiation_isotropic, radiation, V, I, Z = sol

        # @show abs2.(waves["o2@0+"]) ./ flux_radiation["o2+"]
        # @show abs2.(waves["o2@0+"]) ./ flux["o2"]

        fields = OrderedDict(string.(eachindex(fields)) .=> fields)

        if ordering == :frequency
            fields = reverse.(fields)
            waves = vmap(reverse, waves)
            flux = vmap(reverse, flux)
            V = vmap(reverse, V)
            I = vmap(reverse, I)
            Z = vmap(reverse, Z)
            radiation_isotropic = vmap(a -> reverse(a, dims=ndims(a)), radiation_isotropic)
            radiation = vmap(a -> reverse(a, dims=ndims(a)), radiation)
        end
        dBi = vmap(radiation_isotropic) do x
            10log10.(max.(x, 1f-30))
        end
        global solution = (; waves, flux, fields, dBi, radiation, radiation_isotropic, V, I, Z) |> cpu
        writejsonnp(SOL, solution)
    elseif study_type == "inverse_design"
        allowscalar(true)
        println("starting optimization... first iter will be slow due to adjoint compilation.")
        println("")

        η = 1.
        stop = false
        ls = []
        A = sum(prod.(size.(models)))
        g0 = nothing

        momentum0 = momentum
        CONTRAST = linear_interpolation([lowloss, 1] |> FF, [1, contrast] |> FF, extrapolation_bc=Flat())

        ENV["AD"] = 0
        sol = foo(models; prob)
        l = loss(sol; prob)
        @show prob[:umax] = FF(1.1) * sol[:umax]
        ENV["AD"] = 1
        @debug (; lowloss, stoploss)
        println()

        for iter = 1:iters+1
            isckpt = iter == 1 || iter % ckptat == 0 || stop

            if stop
                l = loss(models; prob)
            else
                println("[$iter]")

                momentum = (iter - 1) / iter * momentum0 |> FF
                contrast = max(contrast, CONTRAST(l))
                [(m.contrast = contrast) for m = models]
                @debug (; contrast, momentum)

                @time "gradient" l, (g,) = Zygote.withgradient(models) do ms
                    loss(ms; prob)
                end

                @assert !isnothing(g)
                @assert !any(isnan, g[1].p)
                @assert maximum(abs.(g[1].p)) > 0
            end

            push!(ls, l)
            @debug (; ls)
            stop = l < stoploss || (iter == iters)


            if !stop
                c = 0.1
                maxdA = max(0.1c, c * l)
                if iter == 1
                    mindA = 0.5maxdA
                else
                    mindA = 0.00
                    if ls[end] < ls[end-1]
                        η *= 1.1
                    else
                        η /= 1.15
                    end
                    η = FF(η)

                    g = momentum * g0 + (1 - momentum) * g
                end

                @debug (; mindA, maxdA, η)
                undershot = overshot = false
                models0 = deepcopy(models)
                As0 = [m() for m = models]
                dA = As = nothing
                j = 0
                r = FF(1.05)
                while 0 < η < 0.9floatmax(FF) && j < 2000
                    dx = -η * g
                    dxmax = maximum(abs.(dx[1].p))
                    for (m, m0, dm) = zip(models, models0, dx)
                        m.p .= dm.p + m0.p
                    end
                    As = [m() for m = models]
                    dA = sum(sum.(abs, As - As0)) / A
                    @debug (; dA, η, dxmax)
                    if dA < mindA && !overshot
                        undershot = true
                        η *= r
                    elseif dA > maxdA && !undershot
                        overshot = true
                        η /= r
                    else
                        println(j)
                        break
                    end
                    j += 1
                end

                g0 = g
                println("updated parameters")
                println("area change fraction: $dA")
            end

            if isckpt
                println("saving checkpoint...")
                makemovie(path)
                CKPT = joinpath(path, "checkpoints", replace(string(now()), ':' => '_', '.' => '_'))
                mkpath(CKPT)
                for (_d, d) = zip(designs, prob(:design_instances))
                    for p = (path, CKPT)
                        m = d.model
                        p = joinpath(p, "designs")
                        mkpath(p)
                        npzwrite(joinpath(p, "params_$(d.name).npy"), m.p)
                        a = m()
                        npzwrite(joinpath(p, "design_$(d.name).npy"), resize(a, round.(Int, d.L * λ / pixel_size)))
                    end
                end
                # if gpu_backend == "CUDA"
                println()
                CUDA.memory_status()
                # end
            end
            stop && break
            GC.gc()
            println("----\n")
        end

    end
    0
end

