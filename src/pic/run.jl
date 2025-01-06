function picrun(path; array=Array, kw...)
    Random.seed!(1)
    ENV["autodiff"] = "0"
    println("setting up simulation...")
    global PROB_PATH = joinpath(path, "problem.json")
    SOL_PATH = joinpath(path, "solution.json")
    temp = joinpath(path, "temp")

    io = open(PROB_PATH)
    s = read(io, String)
    global prob = JSON.parse(s; dicttype=OrderedDict)
    # merge!(prob, kw)
    for (k, v) = pairs(kw)
        prob[string(k)] = v
    end
    @unpack name, N, approx_2D_mode, dtype, wl, xmargin, ymargin, runs, ports, dl, xs, ys, zs, components, study, zmode, hmode, zmin, zcenter, gpu_backend, magic, framerate, layer_stack, matprops, L, Ttrans, Tss = prob
    if study == "inverse_design"
        @unpack lsolid, lvoid, designs, targets, weights, eta, iters, restart, save_memory, design_config, stoploss = prob
    end
    F = Float32
    alg = :spectral
    alg = nothing
    dtype = lowercase(dtype)
    if contains(dtype, "16") && contains(dtype, "bf")
        F = BFloat16
        println("BFloat16 selected. make sure your GPU supports it. otherwise will be emulated and very slow.")

    elseif contains(dtype, "16")
        F = Float16
        println("Float16 selected. make sure your cpu or GPU supports it. otherwise will be emulated and very slow.")
    end
    λ = wl
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
        a = stack(map(sort(collect(readdir(joinpath(temp, string(k)), join=true)))) do file
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

        ϵ = matprops(material).ϵ
        ϵmin = min(ϵ, ϵmin)
        ϵ3[I...] .*= 1 - a
        ϵ3[I...] .+= a .* ϵ
    end
    ϵ3 = max.(ϵmin, ϵ3)
    ϵ2 = ϵ3[:, :, 1+round((zcenter - zmin) / dl)]

    global models = nothing
    lb = components.device.bbox[1]
    if N == 2
        ϵ = ϵ2
    else
        ϵ = ϵ3
    end
    if study == "inverse_design"
        targets = fmap(F, targets)
        # targets = sortkeys(targets)
        if isfile(SOL_PATH)
            sol = open(SOL_PATH, "r") do f
                JSON.parse(f)
            end
        else
            sol = nothing
        end
        models = [
            begin
                @unpack bbox = design
                L = bbox[2] - bbox[1]
                szd = Tuple(round.(Int, L / dl)) # design region size
                symmetries = [Int(s) + 1 for s = design.symmetries]

                frame = ϵ2
                frame = frame .>= 0.99maximum(frame)
                # frame = nothing
                start = round((bbox[1] - lb) / dl + 1)
                b = Blob(szd; solid_frac=0.95, lsolid=lsolid / dl, lvoid=lvoid / dl, symmetries, F, frame, start, morph=true)
                display(heatmap(b.frame))

                if !isnothing(sol) && !restart
                    println("loading saved design...")
                    b.p .= sol.params[i]
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
                dimsperm = getdimsperm(L3)
                λmodenums = SortedDict([(F(_λ) / λ) => v for (_λ, v) in pairs(wavelength_mode_numbers)])
                push!(sources, Source(center, -L / 2, L / 2, dimsperm, N, approx_2D_mode, center3, -L3 / 2, L3 / 2; λmodenums, label="s$(string(port)[2:end])"))
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

            Monitor(center, -L / 2, L / 2, dimsperm, N, approx_2D_mode, center3, -L3 / 2, L3 / 2; λmodenums, label=port)
        end for (port, m) = SortedDict(run.monitors) |> pairs] for run in runs]

    global run_probs =
        [
            begin
                setup(dl / λ, boundaries, sources, monitors, deltas[1:N] / λ, mode_deltas[1:N-1] / λ, ; approx_2D_mode, array,
                    F, ϵ, ϵ3, deltas3=deltas / λ, λ, temp, Ttrans, Tss)
            end for (i, (run, sources, monitors)) in enumerate(zip(runs, runs_sources, runs_monitors))
        ]

    # error("not implemented")
    t0 = time()
    lb3 = (lb..., zmin)
    # error("not implemented")
    println("compiling simulation code...")
    if study == "sparams"
        @unpack S, sols = calc_sparams(runs, run_probs, lb, dl;
            F, verbose=true, framerate, path)
        plotsols(sols, run_probs, path)
        sol = (; sparam_family(S)...,
            path, study)
        open(SOL_PATH, "w") do f
            write(f, json(cpu(sol)))
        end
    elseif study == "inverse_design"
        ENV["autodiff"] = "1"
        if length(lb) == 3
            if magic != "summersale"
                error("3D inverse design feature must be requested from Luminescent AI info@luminescentai.com")
            end
        end
        model = models[1]
        minchange = 0.001
        maxchange = max(5minchange, 1.2jump(model) / prod(size(model)))
        global opt = AreaChangeOptimiser(model;
            minchange,
            maxchange,
            # opt=Adam(1, (0.8, 0.9)), 
            opt=Momentum(1, 0.5),
        )
        opt_state = Flux.setup(opt, model)
        # error("not implemented")
        println("starting optimization... first iter will be slow due to adjoint compilation.")
        img = nothing
        println("")
        for i = 1:iters
            println("====($i)  ")
            stop = i == iters
            if :phase_shifter == first(keys(targets))
                @time l, (dldm,) = Flux.withgradient(model) do model
                    models = [model]
                    @unpack S, sols = make_pic_sim_prob(runs, run_probs, lb, dl,
                        designs, design_config, models;
                        F, img, alg)#(1)(1)
                    k = keys(S) |> first
                    s = S[k][Symbol("o2@0,o1@0")]

                    @unpack S = make_pic_sim_prob(runs, run_probs, lb, dl,
                        designs, design_config, models;
                        F, img, alg, save_memory, perturb=:ϵ,)#(1)(1)
                    s_ = S[k][Symbol("o2@0,o1@0")]

                    T = abs2(s)
                    dϕ = angle(s_ / s)
                    println("T: $T, dϕ: $dϕ")
                    (exp(T - 1) * dϕ / π)
                    # T * dϕ / π
                end
            else
                function f(model)
                    models = [model]
                    res = calc_sparams(runs, run_probs, lb, dl,
                        designs, design_config, models, ;
                        F, img, alg, save_memory, matprops)
                    # return res
                    @unpack S, sols = res
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
                        _l = sum(err.(yhat, y) / Z,) do x
                            abs(x)
                        end
                        println("$(k) loss: $_l ")
                        l += _l
                    end
                    # println("$l\n")
                    println()
                    l
                end

                @time global l, (dldm,) = Flux.withgradient(f, model)

            end
            @assert !isnothing(dldm)
            if !isnothing(stoploss) && l < stoploss
                println("\nLoss below threshold, stopping optimization.")
                stop = true
            end
            if true# i == 1 || i % 2 == 0 || stop
                println("saving checkpoint...")
                ckptpath = joinpath(path, "checkpoints", replace(string(now()), ':' => '_', '.' => '_'))

                mkpath(ckptpath)
                models = [model]
                for (i, (m, design)) = enumerate(zip(models, designs))
                    # a = Gray.(m() .< 0.5)

                    # Images.save(joinpath(ckptpath, "optimized_design_region_$i.png"), a)
                    # Images.save(joinpath(path, "optimized_design_region_$i.png"), a)
                end
                plotsols(sols, run_probs, (path, ckptpath))

                sol = (;
                    sparam_family(S)...,
                    optimized_designs=[m() .> 0.5 for m in models], dl,
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
                println(json(sol.T, 4))
            end
            if stop
                break
            end
            opt.minchange = minchange * (1 + 1l)
            opt.maxchange = maxchange * (1 + 4l)
            Jello.update_loss!(opt, l)
            Flux.update!(opt_state, model, dldm)# |> gpu)
            GC.gc()
            println("====\n")

        end
        if framerate > 0
            make_pic_sim_prob(runs, run_probs, lb, dl,
                designs, design_config, models;
                F, img, alg, framerate, path)
        end
        println("Done in $(time() - t0) .")

    end
    if length(sol.T) > 1
        println("wavelengths may have been adjusted to facilitate simulation.")
    end
    println("T-params: ")
    println(JSON.json(sol.T, 4))
    # println(sol)
    sol
end