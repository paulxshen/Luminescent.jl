using Optimisers
using SparseArrays
using Flux: Adam
for T in (:Float32, :Float16, :Float64)
    @eval Optimisers.maywrite(::CUDA.CUSPARSE.CuSparseMatrixCSC{$T,Int32}) = true
    @eval Optimisers.maywrite(::SparseArrays.SparseMatrixCSC{$T,Int32}) = true
end

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

function write_sparams(runs, run_probs, origin, Δ,
    designs=nothing, design_config=nothing, models=nothing;
    alg=nothing, save_memory=false, verbose=false, perturb=nothing, framerate=0, path="", kw...)
    F = run_probs[1].F
    dx = Δ[1]
    sols = [
        begin
            prob[:_geometry] = make_geometry(models, origin, dx, prob._geometry, designs, design_config; F, perturb)
            #@debug typeof(prob.u0.E.Ex), typeof(prob.geometry.ϵ)
            sol = solve(prob; alg, save_memory, verbose, framerate, path)
        end for (i, prob) in enumerate(run_probs)
        # end for (i, prob) in enumerate(run_probs)
    ]
    # S = sols[1]("a+", 1) |> abs2
    # return (; S, sols)

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
                    s = "$monitor_port@$mn," * "$source_port@$source_mn"
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
    isnothing(models) && return geometry
    ratio = models[1].ratio
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

            for (m, design) in zip(models, designs)
                mask = m((x, r) -> x, v, f)
                # mask = m() * (f - v) + v

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

function plotsols(sols, probs, path)
    for (i, (prob, sol)) in enumerate(zip(probs, sols))
        # try
        @unpack u, p, _p = sol |> cpu
        @unpack monitor_instances, source_instances, Δ, dl, λ, ratio = prob |> cpu
        u = u.Hz
        g = _p.ϵ
        u = imresize(u, Tuple(round.(Int, size(u) .* ratio)))
        g = imresize(g, size(u))

        plt = quickie(u, g; dx=dl, λ, monitor_instances, ratio, source_instances,)
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

    global prob = load(PROB_PATH)
    @load PROB_PATH name N dtype xmargin ymargin zmargin dx0 source_margin port_source_offset source_portsides nonsource_portsides runs ports dl dx dy dz components study mode_solutions eps_2D eps_3D hmode zmin thickness zcore zcenter gpu_backend magic framerate
    F = Float32
    alg = :spectral
    alg = nothing
    if contains(dtype, "16")
        F = Float16
        println("Float16 selected. make sure your cpu or GPU supports it. otherwise will be emulated and very slow.")
        # println("Float16 not supported yet, will be in future release.")
    end
    λc = median(load(PROB_PATH)[:wavelengths])

    global eps_2D = convert.(F, stack(stack.(eps_2D)))'
    global eps_3D = convert.(F, permutedims(stack(stack.(eps_3D)), (3, 2, 1)))


    Δ = if N == 2
        [dx, dy]
    else
        [dx, dy, dz]
    end
    ratio = Int.(Δ / dl)

    models = nothing
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
        # for (mode, mode1) in zip(ms.modes, ms.modes1)
        for (mode) in ms.modes
            global _mode = mode
            mode = NamedTuple([k => complex.(stack.(F(v))...) for (k, v) in mode |> pairs])
            # mode1 = (; [k => complex.([convert.(F, v) for v = v]...) for (k, v) in mode1 |> pairs]...)
            if N == 2
                # mode = mode1
                mode = collapse_mode(mode,)
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
    #  a = mode_solutions
    # modal source
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
                        push!(sources, ModalSource(t -> cispi(2t * λc / λ), mode, c, n, tangent, L; meta=(; port)))

                    end
                end
                # sig[:calibrated_modes] = []
            end
            sources
        end for run in runs
    ]

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
            ModalMonitor(wavelength_modes, c, n, tangent, L; meta=(; port))
        end
        for (port, m) = run.monitors |> pairs] for run in runs]

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
        run_probs = gpu.(run_probs)
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
            write_sparams(runs, run_probs, origin, Δ,
                designs, design_config, models;
                F, img, alg, framerate, path)
        end
        println("Done in $(time() - t0) .")

    end
    sol
end