
function picrun(path; kw...)
    Random.seed!(1)
    println("setting up simulation...")
    PROB_PATH = joinpath(path, "problem.bson")
    SOL_PATH = joinpath(path, "solution.json")

    verbose = false

    global prob = load(PROB_PATH)
    @load PROB_PATH name N dtype xmargin ymargin zmargin dx0 source_margin port_source_offset source_portsides nonsource_portsides runs ports dl dx dy dz components study mode_solutions hmode zmin thickness zcore zcenter gpu_backend magic framerate layer_stack materials L
    for (k, v) in pairs(kw)
        @show k, v
        @eval $k = $v
    end

    F = Float32
    alg = :spectral
    alg = nothing
    if contains(dtype, "16")
        F = Float16
        println("Float16 selected. make sure your cpu or GPU supports it. otherwise will be emulated and very slow.")
        # println("Float16 not supported yet, will be in future release.")
    end
    λc = median(load(PROB_PATH)[:wavelengths])

    global eps_2D = nothing
    global sz = Int.(L .÷ dl)
    global eps_3D = zeros(F, Tuple(sz))
    global layer_stack
    layer_stack = sort(collect(pairs(layer_stack)), by=kv -> -kv[2][:mesh_order]) |> OrderedDict
    global ϵmin = Inf
    for (k, v) = pairs(layer_stack)
        a = stack(map(sort(collect(readdir(joinpath(path, "temp", string(k)), join=true)))) do file
            a = F.(Gray.(FileIO.load(file)))
            reverse(a', dims=2)
        end)
        @unpack material, thickness = v

        start = 1 + floor.([(v.origin / dl)..., (v.zmin - zmin) / dl])
        I = range.(start, start + size(a) - 1)
        ϵ = materials[Symbol(material)].epsilon
        ϵmin = min(ϵ, ϵmin)
        eps_3D[I...] .*= 1 - a
        eps_3D[I...] .+= a .* ϵ
    end
    eps_3D = max.(ϵmin, eps_3D)
    eps_2D = eps_3D[:, :, round(Int, (zcenter - zmin) / dl)]
    heatmap(eps_2D) |> display

    error()
    # s = run_probs[1].source_instances[1]
    # ab =Functors.functor(s)
    # a = gpu(s)
    # aa = gpu(s.g.Jy)
    Δ = if N == 2
        [dx, dy]
    else
        [dx, dy, dz]
    end
    ratio = round.(Int, Δ / dl)

    global models = nothing
    # heatmap(eps_2D) |> display
    # GLMakie.volume(eps_3D) |> display
    polarization = :TE
    global origin = components.device.bbox[1]
    if N == 2
        ϵ = eps_2D
    else
        ϵ = eps_3D
    end
    if study == "inverse_design"
        @load PROB_PATH designs targets weights eta iters restart save_memory design_config stoploss
        targets = fmap(F, targets)
        prob = load(PROB_PATH)
        if isfile(SOL_PATH)
            sol = open(SOL_PATH, "r") do f
                JSON.parse(f)
            end
        else
            sol = nothing
        end
        models = [
            begin
                @unpack init, bbox = d
                L = bbox[2] - bbox[1]
                szd = Tuple(round.(Int, L / dx)) # design region size
                symmetries = [length(string(s)) == 1 ? Int(s) + 1 : s for s = d.symmetries]
                o = round.(Int, (bbox[1] - origin) / dx) + 1

                lvoid = d.lvoid / dx
                lsolid = d.lsolid / dx
                margin = maximum(round.(Int, (lvoid, lsolid)))
                ϵ2 = downsample(eps_2D, ratio[1])
                frame = ϵ2[range.(o - margin, o + szd + margin - 1)...]
                frame = frame .== maximum(frame)
                # display(heatmap(frame))
                b = Blob(szd;
                    init, lvoid, lsolid, symmetries, F, frame, ratio=ratio[1])

                if !isnothing(sol) && !restart
                    println("loading saved design...")
                    b.a .= sol.params[i] |> typeof(b.a)
                end
                b
            end for (i, d) = enumerate(designs)
        ]
    end

    boundaries = [] # unspecified boundaries default to PML

    device = 0
    device, dx, λc =
        convert.(F, (device, dx, λc))
    # guess = convert.(F,guess)
    global mode_solutions
    for ms = mode_solutions
        if !haskey(ms, :calibrated_modes)
            ms[:_modes] = []
            ms[:calibrated_modes] = []
        end
        for (mode, mode1D) in zip(ms.modes, ms.modes1D)
            # for (mode) in ms.modes
            global _mode = mode
            mode = NamedTuple([k => complex.(stack.(F(v))...) for (k, v) in mode |> pairs])
            mode1D = (; [k => complex.([convert.(F, v) for v = v]...) for (k, v) in mode1D |> pairs]...)

            if N == 2
                mode = mode1D
                # mode = collapse_mode(mode,)
            end
            push!(ms[:_modes], mode)
            mode = kmap(mode) do a
                # imresize(a, Tuple(round.(Int, size(a) * dl ./ (N == 2 ? Δ[1] : Δ[2:3]))))
                downsample(a, N == 2 ? ratio[1] : ratio[2:3])
            end
            mode = keepxy(mode)
            push!(ms[:calibrated_modes], mode)
        end
    end

    global runs_sources = [
        begin
            sources = []
            for (port, sig) = run.sources |> pairs
                @unpack center, wavelength_mode_numbers = sig
                for λ in keys(wavelength_mode_numbers)
                    for mn = wavelength_mode_numbers[λ]
                        λ = convert.(F, parse(Float32, string(λ)))
                        i = findfirst(mode_solutions) do v
                            abs(λ - v.wavelength) < 0.001 && string(port) in string.(v.ports)
                        end
                        ms = mode_solutions[i]
                        mode = ms.calibrated_modes[mn+1]
                        # sum(abs.(aa.source_instances[1].g.Jy))
                        # heatmap(abs.(bb.source_instances[1]._g.Jy))
                        n = -sig.normal
                        tangent = [-n[2], n[1]]
                        c = (sig.center - origin) / λc
                        L = [sig.mode_width] / λc
                        if N == 3
                            L = [L..., hmode / λc]
                            n, tangent, = vcat.((n, tangent,), ([0],))
                            c = [c..., (zcenter - zmin) / λc]
                        end
                        push!(sources, ModalSource(t -> cispi(2t * λc / λ), mode, c, n, tangent, L; label="s$(string(port)[2:end])"))
                    end
                end
            end
            sources
        end for run in runs
    ]
    # sort!(runs_sources, by=x -> x.label)

    global runs_monitors = [[
        begin
            c = (m.center - origin) / λc
            n = m.normal
            tangent = [-n[2], n[1]]
            L = [m.mode_width] / λc

            if N == 3
                L = [L..., hmode / λc]
                n, tangent, = vcat.((n, tangent,), ([0],))
            end

            wavelength_modes = SortedDict([
                begin
                    λ = parse(Float32, string(λ))
                    λ = F(λ)
                    λ / λc => [
                        begin
                            i = findfirst(mode_solutions) do v
                                abs(λ - v.wavelength) < 0.001 && string(port) in v.ports && mn < length(v.modes)
                            end
                            ms = mode_solutions[i]
                            mode = ms.calibrated_modes[mn+1]
                        end for mn = 0:maximum(mns)
                    ]
                end
                for (λ, mns) in pairs(m.wavelength_mode_numbers)
            ])
            if N == 3
                c = [c..., (zcenter - zmin) / λc]
            end
            ModalMonitor(wavelength_modes, c, n, tangent, L; label=port)
        end for (port, m) = run.monitors |> pairs] for run in runs]
    # sort!(runs_monitors, by=x -> x.label)


    global run_probs =
        [
            begin

                setup(boundaries, sources, monitors, Δ / λc, dl / λc;
                    F, ϵ, ratio,
                    verbose, λ=λc)
            end for (i, (run, sources, monitors)) in enumerate(zip(runs, runs_sources, runs_monitors))
        ]

    if !isempty(gpu_backend)
        println("using $gpu_backend backend.")
        # Flux.gpu_backend!(gpu_backend)
        if gpu_backend == "CUDA"
            # study == "inverse_design" && CUDA.allowscalar(true)
            @assert CUDA.functional()
            # elseif gpu_backend == "AMDGPU"
            #     using AMDGPU
            # elseif gpu_backend == "Metal"
            #     using Metal
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
    global sols = 0
    t0 = time()
    if study == "sparams"
        println("Computing s-parameters...")
        @unpack S, sols = write_sparams(runs, run_probs, origin, Δ;
            F, verbose=true, framerate, path)
        plotsols(sols, run_probs, path, origin)
        sol = (; sparam_family(S)...,
            path, dx, study)
        open(SOL_PATH, "w") do f
            write(f, json(cpu(sol)))
        end
        println("Done in $(time() - t0) s")
    elseif study == "inverse_design"
        if length(origin) == 3
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
            global sols
            println("($i)  ")
            stop = i == iters
            if :phase_shifter == first(keys(targets))
                @time l, (dldm,) = Flux.withgradient(models) do m
                    @unpack S, sols = write_sparams(runs, run_probs, origin, Δ,
                        designs, design_config, models;
                        F, img, alg)#(1)(1)
                    k = keys(S) |> first
                    s = S[k][Symbol("o2@0,o1@0")]

                    @unpack S = write_sparams(runs, run_probs, origin, Δ,
                        designs, design_config, m;
                        F, img, alg, save_memory, perturb=:ϵ,)#(1)(1)
                    s_ = S[k][Symbol("o2@0,o1@0")]

                    T = abs2(s)
                    dϕ = angle(s_ / s)
                    println("T: $T, dϕ: $dϕ")
                    (exp(T - 1) * dϕ / π)
                    # T * dϕ / π
                end
            else
                @time l, (dldm,) = Flux.withgradient(models) do models
                    # sols = get_sols(runs, run_probs,  path, origin, Δ,
                    @unpack S, sols = write_sparams(runs, run_probs, origin, Δ,
                        designs, design_config, models;
                        F, img, alg, save_memory)
                    l = 0
                    for k = keys(targets)
                        y = targets[k]
                        err = -
                        if :phasediff == k
                            ŷ = namedtuple([
                                λ => namedtuple([
                                    ps =>
                                        begin
                                            ks = ignore_derivatives() do
                                                ps = split(string(ps), ",")
                                                ks = keys(S[λ])
                                                [ks[findfirst(k -> startswith(string(k), p), ks)] for p = ps]
                                            end
                                            s1, s2 = [S[λ][k] for k in ks]
                                            angle(s1 / s2)
                                        end
                                    for ps = keys(targets[k][λ])])
                                for λ = keys(targets[k])])
                            err = (x, y) -> angle(cis(x) / cis(y))
                            ŷ = flatten(ŷ)
                            y = flatten(y)
                            Z = length(y) * convert.(F, 2π)
                        else
                            if :sparams == k
                                # S = get_sparams(sols)
                                ŷ = S
                            elseif :tparams == k
                                ŷ = fmap(abs2, S)
                            end

                            # global a1 = ŷ
                            # global a2 = y
                            # println("ŷ: $ŷ")
                            # println("y: $y")

                            ŷ = [[ŷ(λ)(k) for k = keys(y[λ])] for λ = keys(y)]
                            ŷ = flatten(ŷ)
                            y = flatten(y)
                            Z = sum(abs, y)
                        end
                        _l = sum(abs, err.(ŷ, y),) * weights(k) / Z
                        println("$(k) loss: $_l ")
                        # ignore_derivatives() do
                        #     println(json(ŷ))
                        # end
                        l += _l
                    end
                    println("    weighted total loss $l")
                    l
                end
            end

            # if i == 1
            #     best0 = best = l
            # end

            # if l < best
            #     best = l
            # end
            if !isnothing(stoploss) && l < stoploss
                println("Loss below threshold, stopping optimization.")
                stop = true
            end
            # if i % 15 == 0
            #     if best - best0 > -0.01
            #         println("Loss stagnating, stopping optimization.")
            #         stop = true
            #     else
            #         best0 = best
            #     end
            # end
            println("")

            if i == 1 || i % 2 == 0 || stop
                println("saving checkpoint...")
                ckptpath = joinpath(path, "checkpoints", replace(string(now()), ':' => '_', '.' => '_'))

                mkpath(ckptpath)
                for (i, (m, d)) = enumerate(zip(models, designs))
                    # a = Gray.(m() .< 0.5)

                    # Images.save(joinpath(ckptpath, "optimized_design_region_$i.png"), a)
                    # Images.save(joinpath(path, "optimized_design_region_$i.png"), a)
                end
                plotsols(sols, run_probs, (path, ckptpath), origin)

                sol = (;
                    sparam_family(S)...,
                    # optimized_designs=[m() .> 0.5 for m in models],
                    # params=getfield.(models, :a),
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
            # if stop
            #     break
            # end
            global aaa = dldm
            Flux.update!(opt_state, models, dldm)# |> gpu)
        end
        if framerate > 0
            write_sparams(runs, run_probs, origin, Δ,
                designs, design_config, models;
                F, img, alg, framerate, path)
        end
        println("Done in $(time() - t0) .")

    end
    sol
end