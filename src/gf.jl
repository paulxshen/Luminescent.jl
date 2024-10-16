using Optimisers
using SparseArrays
using Flux: Adam
Optimisers.maywrite(::CUDA.CUSPARSE.CuSparseMatrixCSC{Float32,Int32}) = true
Optimisers.maywrite(::CUDA.CUSPARSE.CuSparseMatrixCSC{Float16,Int32}) = true
Optimisers.maywrite(::SparseArrays.SparseMatrixCSC{Float32,Int32}) = true
Optimisers.maywrite(::SparseArrays.SparseMatrixCSC{Float16,Int32}) = true
Optimisers.maywrite(::SparseArrays.SparseMatrixCSC{TrackedFloat32,Int32}) = true

function julia_main()::Cint
    if !isempty(ARGS)
        gfrun(ARGS[1])
    end
    return 0
end

function lastrun(; name=nothing, study=nothing, wd="runs")
    path = joinpath(pwd(), wd)
    l = filter(isdir, readdir(path, join=true))
    sort!(l, by=p -> Dates.unix2datetime(mtime(p)), rev=true)

    if !isnothing(study)
        for p = l
            try
                open(joinpath(p, "solution.json")) do f
                    JSON.parse(f)["study"]
                end == study && return p
            catch e
                println(e)
            end
        end
    end
    return l[1]
end

function write_sparams(runs, run_probs, origin, dx,
    designs=nothing, design_config=nothing, models=nothing;
    alg=nothing, save_memory=false, verbose=false, perturb=nothing = nothing, framerate=0, path="", kw...)
    F = run_probs[1].F

    sols = [
        begin
            prob[:geometry] = make_geometry(models, origin, dx, prob.geometry, designs, design_config; F, perturb)
            #@debug typeof(prob.u0.E.Ex), typeof(prob.geometry.ϵ)
            sol = solve(prob; alg, save_memory, verbose, framerate, path)
        end for (i, prob) in enumerate(run_probs)
        # end for (i, prob) in enumerate(run_probs)
    ]
    S = sols[1]("a+", 1) |> abs2
    return (; S, sols)

    ulims = sols[1].ulims
    # return sol
    coeffs = OrderedDict()
    for (sol, run) in zip(sols, runs)
        sources = run.sources |> Porcupine.values
        monitors = run.monitors |> Porcupine.values
        source_port = first(sources).port
        # source_mn = first(sources).wavelength_mode_numbers(1)[1]
        source_mn = first(sources).wavelength_mode_numbers |> Porcupine.first |> Porcupine.first
        for (m, monitor) = enumerate(monitors)
            for (w, λ) = enumerate(keys(monitor.wavelength_mode_numbers))
                for mn = monitor.wavelength_mode_numbers[λ]
                    monitor_port = monitor.port
                    λ = Symbol(λ)
                    if !haskey(coeffs, λ)
                        coeffs[λ] = OrderedDict()
                    end
                    s = "o$monitor_port@$mn," * "o$source_port@$source_mn"
                    s = Symbol(s)
                    coeffs[λ][s] = (sol("a+", m, w, mn), sol("a-", m, w, mn))
                end
            end
        end
    end
    # return coeffs(1)(1)[1] |> abs2

    S = OrderedDict([λ => OrderedDict([k => begin
        s = ignore() do
            split(string(k), ",")[2]
        end
        # Symbol(
        coeffs[λ][k][1] / coeffs[λ][Symbol("$s,$s")][2]
    end for (k) = keys(coeffs[λ])]) for (λ) = keys(coeffs)])
    # if source_mn == mn == 0
    #     coeffs[λ]["$monitor_port,$source_port")] = v
    # end
    return (; S, sols)
end

function make_geometry(models, origin, dx, geometry, designs, design_config; F=Float32, perturb=nothing)
    return (; μ=geometry.μ, σ=geometry.σ, m=geometry.m, ϵ=2.5mean(models[1]) * ones(F, size(geometry.ϵ)))
    isnothing(models) && return geometry
    ratio = 1
    # ratio = models[1].ratio
    dx /= ratio
    namedtuple([k => begin
        if k in keys(design_config.fill)
            f = design_config.fill[k] |> F
            v = design_config.void[k] |> F
            if perturb == k
                f *= convert.(F, 1.001)
            end

            b = Zygote.Buffer(geometry[k])
            copyto!(b, geometry[k])

            for (mask, design) in zip(models, designs)
                # global mask = m((x, r) -> x, v, f)
                mask = mask * (f - v) + v

                o = round.(Int, (design.bbox[1] - origin) / dx) + 1
                if ndims(b) == 3
                    o = [o..., 1 + round(Int, (zcore - zmin) / dx)]
                    mask = stack(fill(mask, round(Int, thickness / dx)))
                end
                # b[range.(o, o .+ size(mask) .- 1)...] = mask
                b[[i:j for (i, j) = zip(o, o .+ size(mask) .- 1)]...] = mask
            end
            copy(b)
        else
            geometry[k]
        end
    end for k = keys(geometry)])
end

function plotsols(sols, probs, path; ratio=1)
    for (i, (prob, sol)) in enumerate(zip(probs, sols))
        # try
        @unpack u, p, = sol |> cpu
        @unpack monitor_instances, source_instances, dx, λ = prob |> cpu
        a = u.Hz
        # g = _p.ϵ

        plt = quickie(a, ; dx, λ, monitor_instances, source_instances,)
        display(plt)

        if !isa(path, Base.AbstractVecOrTuple)
            path = (path,)
        end
        for path = path
            try
                CairoMakie.save(joinpath(path, "run_$i.png"), plt,)
            catch e
                println("save plot failed")
                println(e)
            end
        end
    end
end

# using AbbreviatedStackTraces
# global virgin, stop, best, best0, sparams0
function gfrun(path; kw...)
    Random.seed!(1)
    println("setting up simulation...")
    PROB_PATH = joinpath(path, "problem.bson")
    SOL_PATH = joinpath(path, "solution.json")

    calibrate = true
    model_name = nothing # if load saved model
    tol = 1e-3
    dosave = false
    verbose = false

    prob = load(PROB_PATH)
    @load PROB_PATH name N dtype margin zmargin dx0 source_margin port_source_offset source_portsides nonsource_portsides runs ports dx components study mode_solutions eps_2D eps_3D mode_height zmin thickness zcore gpu_backend magic framerate ratio
    ratio = 1
    F = Float64
    # F = TrackedFloat32

    alg = :spectral
    alg = nothing
    if contains(dtype, "16")
        F = Float16
        println("Float16 selected. make sure your cpu or GPU supports it. otherwise will be emulated and very slow.")
        # println("Float16 not supported yet, will be in future release.")
    end
    λc = median(load(PROB_PATH)[:wavelengths])

    global eps_2D = convert.(F, stack(stack.(eps_2D)))'
    eps_3D = convert.(F, permutedims(stack(stack.(eps_3D)), (3, 2, 1)))
    eps_2D = eps_2D[range.(1, floor.(Int, size(eps_2D) / ratio) * ratio)...]
    eps_3D = eps_3D[range.(1, floor.(Int, size(eps_3D) / ratio) * ratio)...]


    port_source_offset = trim(port_source_offset, dx)
    margin = trim(margin, dx)
    zmargin = trim(zmargin, dx)
    source_margin = trim(source_margin, dx)

    n_source_portsides = round(Int, (port_source_offset + source_margin) / dx)
    n_nonsource_portsides = round(Int, (source_margin) / dx)
    # p = m + n
    # origin = components.device.bbox[1] - dx * p
    # eps_3D = pad(eps_3D, :replicate, [p, p, 0])
    # ϵ2 = pad(ϵ2, :replicate, p)
    models = nothing
    lr = n_source_portsides * source_portsides + n_nonsource_portsides * nonsource_portsides
    eps_2D = pad(eps_2D, :replicate, ratio * lr...)
    if N == 3
        eps_3D = pad(eps_3D, :replicate, ratio * lr...)
    end
    # ϵ2 = downsample(eps_2D, ratio)

    # global invϵ2 = tensorinv(eps_2D, ratio)
    # heatmap(eps_2D) |> display
    if N == 3
        # ϵ3 = downsample(eps_3D, ratio)
        # invϵ3 = tensorinv(eps_3D, ratio)
    end
    # heatmap(ϵ2) |> display
    # GLMakie.volume(eps_3D) |> display
    #  ϵbase ϵclad ϵcore hbase hwg hclad
    ports, = (ports,) .|> SortedDict
    polarization = :TE

    # _dx = dx / λc
    origin = components.device.bbox[1] - dx * lr[1]

    nz = round(Int, zmargin / dx)
    push!(lr[1], 0)
    push!(lr[2], 0)
    if N == 2
        ϵ = eps_2D
        # ϵ = ϵ2
        # invϵ = invϵ2
        # heatmap(ϵ2) |> display
    else
        ϵ = eps_3D
        # ϵ = ϵ3
        # invϵ = invϵ3
    end
    sz = size(ϵ) .÷ ratio
    ϵmax = maximum(ϵ)
    μ = 1


    # heatmap(ϵ3[round(size(ϵ3, 1) / 2), :, :]) |> display
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
                ϵ2 = downsample(eps_2D, ratio)
                frame = ϵ2[range.(o - margin, o + szd + margin - 1)...]
                frame = frame .== maximum(frame)
                # display(heatmap(frame))
                b = Blob(szd;
                    init, lvoid, lsolid, symmetries, F, frame)

                if !isnothing(sol) && !restart
                    println("loading saved design...")
                    b.a .= sol.params[i] |> typeof(b.a)
                end
                b
                ones(F, szd)
            end for (i, d) = enumerate(designs)
        ]
    end

    boundaries = [] # unspecified boundaries default to PML

    device = 0
    device, dx, λc =
        convert.(F, (device, dx, λc))
    # guess = convert.(F,guess)
    for ms = mode_solutions
        if !haskey(ms, :calibrated_modes)
            ms[:calibrated_modes] = []
        end
        for (mode, mode1) in zip(ms.modes, ms.modes1)
            mode = (; [k => complex.([convert.(F, v) for v = stack.(v)]...) for (k, v) in mode |> pairs]...)
            mode1 = (; [k => complex.([convert.(F, v) for v = v]...) for (k, v) in mode1 |> pairs]...)

            if N == 2
                mode = mode1
                # mode = collapse_mode(mode,)
            end
            # mode = fmap(mode) do x
            #     downsample(x, 2)
            # end
            mode = keepxy(mode)
            push!(ms[:calibrated_modes], mode)
        end
    end
    #  a = mode_solutions
    # modal source
    runs_sources = [
        begin
            d = run.d
            sources = []
            _modes = []
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
                        c = (sig.center - origin + port_source_offset * sig.normal) / λc
                        L = [sig.mode_width] / λc
                        if d == 3
                            L = [L..., mode_height / λc]
                            n, tangent, = vcat.((n, tangent,), ([0],))
                            c = [c..., (ms.zcenter - zmin) / λc]
                        end
                        push!(sources, ModalSource(t -> cispi(2t * λc / λ), mode, c, n, tangent, L; meta=(; port)))

                    end
                end
                # sig[:calibrated_modes] = []
            end
            sources
        end for run in runs
    ]

    runs_monitors = [[
        begin
            d = run.d
            port = parse(Int, string(port))
            c = (m.center - origin) / λc
            n = m.normal
            tangent = [-n[2], n[1]]
            L = [m.mode_width] / λc

            if d == 3
                L = [L..., mode_height / λc]
                n, tangent, = vcat.((n, tangent,), ([0],))
            end

            zcenter = nothing
            wavelength_modes = SortedDict([
                begin
                    λ = parse(Float32, string(λ))
                    λ = F(λ)
                    λ / λc => [
                        begin
                            i = findfirst(mode_solutions) do v
                                abs(λ - v.wavelength) < 0.001 && port in v.ports && mn < length(v.modes)
                            end
                            ms = mode_solutions[i]
                            if isnothing(zcenter)
                                zcenter = ms.zcenter
                            end
                            mode = ms.calibrated_modes[mn+1]
                        end for mn = 0:maximum(mns)
                    ]
                end
                for (λ, mns) in pairs(m.wavelength_mode_numbers)
            ])
            if d == 3
                c = [c..., (zcenter - zmin) / λc]
            end
            ModalMonitor(wavelength_modes, c, n, tangent, L; meta=(; port))
        end
        for (port, m) = run.monitors |> pairs] for run in runs]

    run_probs =
        [
            begin

                setup(boundaries, sources, monitors, dx / λc, sz;
                    F, ϵ,
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
        run_probs = gpu.(run_probs)
        models = models |> gpu
    else
        println("using CPU backend.")
    end

    virgin = true
    # error()
    # Base.iterate(s::Symbol) = println(s)

    t0 = time()
    if study == "sparams"
        # @info "Computing s-parameters..."
        println("Computing s-parameters...")
        # sparams = write_sparams(img="", alg=false, verbose=true)
        @unpack S, sols = write_sparams(runs, run_probs, origin, dx;
            F, verbose=true, framerate, path)
        plotsols(sols, run_probs, path;)
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
        # save_memory = true
        S = sparams0 = 0
        # eta *= 20
        eta = 0.01
        iters = 2
        opt = Flux.Adam(eta)
        opt_state = Flux.setup(opt, models)
        println("starting optimization... first iter will be slow due to adjoint compilation.")
        img = nothing
        best = best0 = 0
        # S = write_sparams(runs, run_probs, origin, dx,
        #     designs, design_config, models;
        #     F, img, alg, with=true,)
        # heatmap(_as[3])
        #  ass = gradient(models) do models
        #     write_sparams(runs, run_probs,  path, origin, dx,
        #         designs, design_config, models;
        #         F, img, alg, save_memory)
        # end
        # error()

        # prob = run_probs[1]
        #  aaaaa = gradient(_geometry0) do geometry
        #     # solve(prob, geometry; alg, save_memory, verbose).forward_mode_powers[1][1][1]
        #     solve(prob, geometry; alg, save_memory, verbose)
        # end
        # error()
        println("")
        sols = 0
        for i = 1:iters
            println("($i)  ")
            stop = i == iters
            # global virgin, stop, best, best0, sparams0
            if :phase_shifter == first(keys(targets))
                @time l, (dldm,) = Flux.withgradient(models) do m
                    @unpack S, sols = write_sparams(runs, run_probs, origin, dx,
                        designs, design_config, models;
                        F, img, alg)#(1)(1)
                    k = keys(S) |> first
                    s = S[k][Symbol("o2@0,o1@0")]

                    @unpack S = write_sparams(runs, run_probs, origin, dx,
                        designs, design_config, m;
                        F, img, alg, save_memory, perturb=:ϵ,)#(1)(1)
                    s_ = S[k][Symbol("o2@0,o1@0")]

                    T = abs2(s)
                    dϕ = angle(s_ / s)
                    println("T: $T, dϕ: $dϕ")
                    S = S
                    # (exp(T - 1) * dϕ / π)
                    T * dϕ / π
                end
            else
                @time l, (dldm,) = Flux.withgradient(models) do models
                    # sols = get_sols(runs, run_probs,  path, origin, dx,
                    @unpack S, sols = write_sparams(runs, run_probs, origin, dx,
                        designs, design_config, models;
                        F, img, alg, save_memory)
                    return S
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
                plotsols(sols, run_probs, (path, ckptpath);)

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
            write_sparams(runs, run_probs, origin, dx,
                designs, design_config, models;
                F, img, alg, framerate, path)
        end
        println("Done in $(time() - t0) .")

    end
    sol
end