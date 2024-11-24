function picrun(path; kw...)
    Random.seed!(1)
    println("setting up simulation...")
    global PROB_PATH = joinpath(path, "problem.json")
    SOL_PATH = joinpath(path, "solution.json")

    io = open(PROB_PATH)
    s = read(io, String)
    global prob = JSON.parse(s)
    @unpack name, N, dtype, xmargin, ymargin, dx0, source_margin, runs, ports, dl, xs, ys, zs, components, study, zmode, hmode, zmin, zcenter, gpu_backend, magic, framerate, layer_stack, materials, L = prob
    if study == "inverse_design"
        @unpack designs, targets, weights, eta, iters, restart, save_memory, design_config, stoploss = prob
    end
    # iters = 15
    # for (k, v) in pairs(kw)
    #     @show k, v
    #     @eval $k = $k
    #     @eval $k = $v
    # end
    # @show eta
    # eta = 0.02
    # N = 2
    # heatmap(debug.mask) 

    F = Float32
    alg = :spectral
    alg = nothing
    if contains(dtype, "16")
        F = Float16
        println("Float16 selected. make sure your cpu or GPU supports it. otherwise will be emulated and very slow.")
    end
    λ = median(prob.wavelengths)
    ticks = [
        begin
            round((v - v[1]) / dl)
        end for v = [xs, ys, zs]
    ]
    spacings = diff.(ticks)
    x = spacings[1][1]
    spacings = [x, x, spacings[3]]
    dx = x * dl
    popfirst!.(ticks)
    deltas = spacings * dl


    ϵ2 = nothing
    sz = round.(L / dl)
    ϵ3 = zeros(F, Tuple(sz))
    layer_stack = sort(collect(pairs(layer_stack)), by=kv -> -kv[2].mesh_order) |> OrderedDict
    ϵmin = Inf
    for (k, v) = pairs(layer_stack)
        a = stack(map(sort(collect(readdir(joinpath(path, "temp", string(k)), join=true)))) do file
            a = F.(Gray.(FileIO.load(file)))
            reverse(a', dims=2)
        end)
        # a = a[Base.OneTo.(min.(size(a), sz))...]
        @unpack material, thickness = v

        start = 1 + round.([(v.origin / dl)..., (v.zmin - zmin) / dl])
        I = range.(start, start + size(a) - 1)
        overhang = max.(last.(I) .- sz, 0)
        a = a[Base.oneto.(size(a) - overhang)...]
        I = range.(start, start + size(a) - 1)

        ϵ = materials(material).epsilon
        ϵmin = min(ϵ, ϵmin)
        ϵ3[I...] .*= 1 - a
        ϵ3[I...] .+= a .* ϵ
    end
    ϵ3 = max.(ϵmin, ϵ3)
    ϵ2 = ϵ3[:, :, 1+round((zcenter - zmin) / dl)]
    # heatmap(ϵ2) |> display
    # GLMakie.volume(ϵ3) |> display
    # error("stop")

    # s = run_probs[1].source_instances[1]
    # ab =Functors.functor(s)
    # a = gpu(s)
    # aa = gpu(s.g.Jy)

    global models = nothing
    lb = components.device.bbox[1]
    if N == 2
        ϵ = ϵ2
    else
        ϵ = ϵ3
    end
    if study == "inverse_design"
        targets = fmap(F, targets)
        if isfile(SOL_PATH)
            sol = open(SOL_PATH, "r") do f
                JSON.parse(f)
            end
        else
            sol = nothing
        end
        models = [
            begin
                @unpack init, bbox = design
                L = bbox[2] - bbox[1]
                szd = Tuple(round.(Int, L / dl)) # design region size
                symmetries = [Int(s) + 1 for s = design.symmetries]

                lvoid = design.lvoid / dl
                lsolid = design.lsolid / dl
                frame = ϵ2
                frame = frame .>= 0.99maximum(frame)
                # frame = nothing
                start = round((bbox[1] - lb) / dl + 1)
                b = Blob(szd; init, lvoid, lsolid, symmetries, F, frame, start)
                # display(heatmap(b.frame))

                if !isnothing(sol) && !restart
                    println("loading saved design...")
                    b.a .= sol.params[i] |> typeof(b.a)
                end
                b
            end for (i, design) = enumerate(designs)
        ]
    end

    boundaries = [] # unspecified boundaries default to PML

    device = 0
    λ = F(λ)
    # guess = convert.(F,guess)
    i = int(v2i(zmode - zmin, deltas[3]))
    j = int(v2i(zmode + hmode - zmin, deltas[3]))
    mode_spacings = [spacings[1][1], adddims(spacings[3][i+1:j], dims=1)]
    mode_deltas = mode_spacings * dl
    global runs = [SortedDict([k => isa(v, AbstractDict) ? SortedDict(v) : v for (k, v) = pairs(run)]) for run in runs]
    global runs_sources = [
        begin
            sources = []
            for (port, sig) = SortedDict(run.sources) |> pairs
                @unpack center, wavelength_mode_numbers = sig
                # i = findfirst(mode_solutions) do v
                #     abs(_λ - v.wavelength) < 0.001 && string(port) in string.(v.ports)
                # end
                # sum(abs.(aa.source_instances[1].g.Jy))
                # heatmap(abs.(bb.source_instances[1]._g.Jy))
                n = -sig.normal
                tangent = [-n[2], n[1]]
                center = (sig.center - lb) / λ
                L = tangent * sig.mode_width / λ
                L3 = [L..., hmode / λ]
                center3 = [center..., (zcenter - zmin) / λ]
                if N == 3
                    L = L3
                    center = center3
                end
                dimsperm = getdimsperm([L..., 1])
                λmodenums = SortedDict([(F(_λ) / λ) => v for (_λ, v) in pairs(wavelength_mode_numbers)])
                push!(sources, Source(center, -L / 2, L / 2, dimsperm, N, center3, -L3 / 2, L3 / 2; λmodenums, label="s$(string(port)[2:end])"))
            end
            sources
        end for run in runs
    ]
    # sort!(runs_sources, by=x -> x.label)

    global runs_monitors = [[
        begin
            center = (m.center - lb) / λ
            n = m.normal
            tangent = [-n[2], n[1]]
            # n, tangent, = vcat.((n, tangent,), ([0],))

            λmodenums = SortedDict([(F(_λ) / λ) => v for (_λ, v) in pairs(m.wavelength_mode_numbers)])

            L = tangent * m.mode_width / λ
            L3 = [L..., hmode / λ]
            center3 = [center..., (zcenter - zmin) / λ]
            if N == 3
                L = L3
                center = center3
            end
            dimsperm = getdimsperm(L3)

            Monitor(center, -L / 2, L / 2, dimsperm, N, center3, -L3 / 2, L3 / 2; λmodenums, label=port)
        end for (port, m) = SortedDict(run.monitors) |> pairs] for run in runs]

    global run_probs =
        [
            begin
                setup(dl / λ, boundaries, sources, monitors, deltas[1:N] / λ, mode_deltas[1:N-1] / λ, ;
                    F, ϵ, ϵ3, deltas3=deltas / λ, λ)
            end for (i, (run, sources, monitors)) in enumerate(zip(runs, runs_sources, runs_monitors))
        ]

    if !isempty(gpu_backend)
        println("using $gpu_backend backend.")
        # Flux.gpu_backend!(gpu_backend)
        if gpu_backend == "CUDA"
            @assert CUDA.functional()
        end
        for prob = run_probs
            for k = keys(prob)
                if k in (:u0, :_geometry, :geometry, :source_instances, :monitor_instances)
                    prob[k] = prob[k] |> gpu
                end
            end
        end
        models = models |> gpu
    else
        println("using CPU backend.")
    end
    t0 = time()
    lb3 = (lb..., zmin)
    if study == "sparams"
        println("Computing s-parameters...")
        @unpack S, sols = write_sparams(runs, run_probs, lb, dl;
            F, verbose=true, framerate, path)
        plotsols(sols, run_probs, path)
        sol = (; sparam_family(S)...,
            path, study)
        open(SOL_PATH, "w") do f
            write(f, json(cpu(sol)))
        end
        println("Done in $(time() - t0) s")
    elseif study == "inverse_design"
        if length(lb) == 3
            if magic != "summersale"
                error("3D inverse design feature must be requested from Luminescent AI info@luminescentai.com")
            end
        end
        opt = Flux.Adam(eta)
        opt_state = Flux.setup(opt, models)
        println("starting optimization... first iter will be slow due to adjoint compilation.")
        img = nothing
        best = best0 = 0
        println("")
        for i = 1:iters
            println("($i)  ")
            stop = i == iters
            if :phase_shifter == first(keys(targets))
                @time l, (dldm,) = Flux.withgradient(models) do models
                    @unpack S, sols, lminloss = write_sparams(runs, run_probs, lb, dl,
                        designs, design_config, models;
                        F, img, alg)#(1)(1)
                    k = keys(S) |> first
                    s = S[k][Symbol("o2@0,o1@0")]

                    @unpack S = write_sparams(runs, run_probs, lb, dl,
                        designs, design_config, models;
                        F, img, alg, save_memory, perturb=:ϵ,)#(1)(1)
                    s_ = S[k][Symbol("o2@0,o1@0")]

                    T = abs2(s)
                    dϕ = angle(s_ / s)
                    println("T: $T, dϕ: $dϕ")
                    (exp(T - 1) * dϕ / π) + lminloss
                    # T * dϕ / π
                end
            else
                @time global l, (dldm,) = Flux.withgradient(models) do models
                    # sols = get_sols(runs, run_probs,  path, lb, deltas,
                    @unpack S, sols = write_sparams(runs, run_probs, lb, dl,
                        designs, design_config, models;
                        F, img, alg, save_memory)
                    l = 0
                    for k = keys(targets)
                        y = targets[k]
                        err = -
                        if :phasediff == k
                            yhat = namedtuple([
                                _λ => namedtuple([
                                    ps =>
                                        begin
                                            ks = ignore_derivatives() do
                                                ps = split(string(ps), ",")
                                                ks = keys(S[_λ])
                                                [ks[findfirst(k -> startswith(string(k), p), ks)] for p = ps]
                                            end
                                            s1, s2 = [S[_λ][k] for k in ks]
                                            angle(s1 / s2)
                                        end
                                    for ps = keys(targets[k][_λ])])
                                for _λ = keys(targets[k])])
                            err = (x, y) -> angle(cis(x) / cis(y))
                            yhat = flatten(yhat)
                            y = flatten(y)
                            Z = length(y) * convert.(F, 2π)
                        else
                            yhat = if "sparams" == string(k)
                                S
                            elseif "tparams" == string(k)
                                fmap(abs2, S)
                            end

                            # global a1 = S, y
                            yhat = [[yhat(_λ)(k) for k = keys(y[_λ])] for _λ = keys(y)]
                            yhat = flatten(yhat)
                            y = flatten(y)
                            Z = sum(abs, y)
                        end
                        _l = sum(abs, err.(yhat, y),) * weights(k) / Z
                        println("$(k) loss: $_l ")
                        l += _l
                    end
                    println("    weighted total loss $l")
                    l
                end
            end

            if !isnothing(stoploss) && l < stoploss
                println("Loss below threshold, stopping optimization.")
                stop = true
            end
            println("")

            if i == 1 || i % 2 == 0 || stop
                println("saving checkpoint...")
                ckptpath = joinpath(path, "checkpoints", replace(string(now()), ':' => '_', '.' => '_'))

                mkpath(ckptpath)
                for (i, (m, design)) = enumerate(zip(models, designs))
                    # a = Gray.(m() .< 0.5)

                    # Images.save(joinpath(ckptpath, "optimized_design_region_$i.png"), a)
                    # Images.save(joinpath(path, "optimized_design_region_$i.png"), a)
                end
                plotsols(sols, run_probs, (path, ckptpath))

                sol = (;
                    sparam_family(S)...,
                    # optimized_designs=[m() .> 0.5 for m in models],
                    params=getfield.(models, :p),
                    designs,
                    design_config, path,
                    dx,
                    study,
                )

                open(joinpath(ckptpath, "solution.json"), "w") do f
                    write(f, json(cpu(sol)))
                end
                open(joinpath(path, "solution.json"), "w") do f
                    write(f, json(cpu(sol)))
                end
            end
            if stop
                break
            end
            Flux.update!(opt_state, models, dldm)# |> gpu)

            # da = Inf
            # α = 1
            # masks0 = [m() for m in models]
            # models0 = deepcopy(models)
            # Flux.update!(opt_state, models, dldm)# |> gpu)
            # models1 = deepcopy(models)
            # while da > 0.5l
            #     for (m, m0, m1) = zip(models, models0, models1)
            #         m.p .= α * m1.p + (1 - α) * m0.p
            #     end
            #     masks = [m() for m in models]
            #     da = sum(masks - masks0) do a
            #         sum(abs, a)
            #     end
            #     da /= sum(masks) do a
            #         prod(size(a))
            #     end
            #     α *= 0.9
            # end
            # @show α
            # println("fractional change in design: $da")
        end
        if framerate > 0
            write_sparams(runs, run_probs, lb, dl,
                designs, design_config, models;
                F, img, alg, framerate, path)
        end
        println("Done in $(time() - t0) .")

    end
    sol
end