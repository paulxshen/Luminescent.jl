using Optimisers
using SparseArrays
using CUDA
Optimisers.maywrite(::CUDA.CUSPARSE.CuSparseMatrixCSC{Float32,Int32}) = true
Optimisers.maywrite(::CUDA.CUSPARSE.CuSparseMatrixCSC{Float16,Int32}) = true
Optimisers.maywrite(::SparseArrays.SparseMatrixCSC{Float32,Int32}) = true
Optimisers.maywrite(::SparseArrays.SparseMatrixCSC{Float16,Int32}) = true

function julia_main()::Cint
    if !isempty(ARGS)

        gfrun(ARGS[1])
    end
    return 0
end
function lastrun(s=nothing, path=joinpath(pwd(), "lumi_runs"))
    l = filter(isdir, readdir(path, join=true))
    sort!(l)
    if isnothing(s)
        return l[end]
    end
    for p = reverse(l)
        try
            open(joinpath(p, "sol.json")) do f
                JSON.parse(f)["study"]
            end == s && return p
        catch e
            println(e)
        end
    end
end

function write_sparams(runs, run_probs, geometry, path, origin, dx,
    designs=nothing, design_config=nothing, models=nothing;
    img=nothing, autodiff=false, compression=false, verbose=false, perturb=nothing, with=false, ls=nothing, kw...)
    F = run_probs[1].F
    geometry = make_geometry(models, origin, dx, geometry, designs, design_config; F, perturb)

    sols = [
        begin
            prob[:geometry] = geometry
            #@debug typeof(prob.u0.E.Ex), typeof(prob.geometry.ϵ)
            sol = solve(prob; autodiff, ls, compression, verbose)

            ignore() do
                if !isnothing(img)
                    if img == ""
                        img = "run$i.png"
                    end
                    try
                        CairoMakie.save(joinpath(path, img), quickie(sol |> cpu, cpu(prob)),)
                    catch e
                        println("plot failed")
                        println(e)
                    end
                end
            end
            sol
        end for (i, prob) in enumerate(run_probs)
        # end for (i, prob) in enumerate(run_probs)
    ]
    # return sols[1].forward_mode_powers[1][1][1]
    # return sols[1].mode_coeffs[1][1][1][1] |> abs2

    ls = sols[1].ls
    # return sol
    coeffs = OrderedDict()
    for (v, run) in zip(sols, runs)
        sources = run.sources |> Porcupine.values
        monitors = run.monitors |> Porcupine.values
        global aaaaaaa = v
        source_port = first(sources).port
        # source_mn = first(sources).wavelength_mode_numbers(1)[1]
        source_mn = first(sources).wavelength_mode_numbers |> Porcupine.first |> Porcupine.first
        for (monitor, v) = zip(monitors, v.mode_coeffs)
            for (λ, v) = zip(monitor.wavelength_mode_numbers |> keys, v)
                for (monitor_mn, v) = zip(monitor.wavelength_mode_numbers[λ], v)
                    monitor_port = monitor.port
                    λ = Symbol(λ)
                    if !haskey(coeffs, λ)
                        coeffs[λ] = OrderedDict()
                    end
                    s = "o$monitor_port@$monitor_mn," * "o$source_port@$source_mn"
                    s = Symbol(s)
                    coeffs[λ][s] = v
                end
            end
        end
    end
    # return coeffs(1)(1)[1] |> abs2

    global sparams = OrderedDict([λ => OrderedDict([k => begin
        s = ignore() do
            split(string(k), ",")[2]
        end
        # Symbol(
        coeffs[λ][k][1] / coeffs[λ][Symbol("$s,$s")][2]
    end for (k) = keys(coeffs[λ])]) for (λ) = keys(coeffs)])
    # if source_mn == monitor_mn == 0
    #     coeffs[λ]["$monitor_port,$source_port")] = v
    # end
    if with
        @show ls
        return sparams, ls
    end
    # return sparams(1)(1) |> abs2
    sparams
end

function make_geometry(models, origin, dx, g, designs, design_config; F=Float32, perturb=nothing)
    isnothing(models) && return g
    # g = deepcopy(g)
    # r = OrderedDict()
    # for k = keys(g)
    #     if haskey(design_config.fill, k)
    #         a = g[k]
    #         for (m, d) in zip(models, designs)
    #             # a = g[k]
    #             mask = m()
    #             f = design_config.fill[k] |> F
    #             v = design_config.void[k] |> F
    #             if perturb == k
    #                 f *= 1.001
    #             end
    #             T = typeof(a)
    #             p = (v .* (1 - mask) + f .* mask) |> F |> T
    #             b = Zygote.Buffer(a)
    #             copyto!(b, a)
    #             xy = round((d.bbox[1] - origin) / dx) + 1
    #             if ndims(a) == 2
    #                 o = xy
    #             else
    #                 o = [xy..., 1 + round((zcore - zmin) / dx)]
    #             end
    #             o = Tuple(o)
    #             b = place!(b, o, p)
    #             a = copy(b)
    #             # g[k] = a
    #         end
    #         r[k] = a
    #     else
    #         r[k] = g[k]
    #     end
    # end
    # r


    # namedtuple(g)
    namedtuple([k => begin
        a = g[k]
        if haskey(design_config.fill, k)
            for (m, d) in zip(models, designs)
                mask = m()
                f = design_config.fill[k] |> F
                v = design_config.void[k] |> F
                if perturb == k
                    f *= 1.001
                end
                T = typeof(a)
                p = (v .* (1 - mask) + f .* mask) |> F
                b = Zygote.Buffer(a)
                copyto!(b, a)
                xy = round((d.bbox[1] - origin) / dx) + 1
                if ndims(a) == 2
                    o = xy
                else
                    o = [xy..., 1 + round((zcore - zmin) / dx)]
                    p = stack(fill(p, round(thickness / dx)))
                end
                o = Tuple(o)
                b = place!(b, o, p)
                a = copy(b)
            end
        end
        a
    end for k = keys(g,)])
end

# using AbbreviatedStackTraces
# using Random
# Random.seed!(0)
# include("main.jl")
# path = lastrun()
# global virgin, stop, best, best0, sparams0
function gfrun(path; kw...)
    println("setting up simulation...")
    PROB_PATH = joinpath(path, "prob.bson")
    SOL_PATH = joinpath(path, "sol.json")

    calibrate = true
    model_name = nothing # if load saved model
    tol = 1e-3
    dosave = false
    verbose = false

    @load PROB_PATH name dtype margin zmargin source_margin Courant port_source_offset portsides runs ports dx components study mode_solutions eps_2D eps_3D mode_height zmin thickness zcore gpu_backend d magic
    F = Float32
    if contains(dtype, "16")
        F = Float16
        println("Float16 selected. make sure your cpu or GPU supports it. otherwise will be emulated and very slow.")
        # println("Float16 not supported yet, will be in future release.")
    end
    λc = median(load(PROB_PATH)[:wavelengths])
    eps_2D = F(stack(stack.(eps_2D)))'
    eps_3D = F(permutedims(stack(stack.(eps_3D)), (3, 2, 1)))
    # heatmap(eps_2D) |> display
    # GLMakie.volume(eps_3D) |> display
    # error()
    #  ϵbase ϵclad ϵcore hbase hwg hclad
    ports, = (ports,) .|> SortedDict
    polarization = :TE

    # _dx = dx / λc
    port_source_offset = whole(port_source_offset, dx)
    margin = whole(margin, dx)
    zmargin = whole(zmargin, dx)
    source_margin = whole(source_margin, dx)

    p = round((port_source_offset + source_margin) / dx)
    np = round(margin / dx)
    # p = m + n
    # origin = components.device.bbox[1] - dx * p
    # eps_3D = pad(eps_3D, :replicate, [p, p, 0])
    # eps_2D = pad(eps_2D, :replicate, p)
    models = nothing
    lr = p * portsides
    eps_2D = pad(eps_2D, :replicate, lr...)

    nz = round(zmargin / dx)
    push!(lr[1], 0)
    push!(lr[2], 0)
    eps_3D = pad(eps_3D, :replicate, lr...)
    # heatmap(eps_3D[round(size(eps_3D, 1) / 2), :, :]) |> display
    origin = components.device.bbox[1] - dx * (p * portsides[1])
    if study == "inverse_design"
        @load PROB_PATH designs targets weights preset eta iters restart compression design_config contrast
        targets = fmap(F, targets)
        prob = load(PROB_PATH)
        minloss = haskey(prob, :minloss) ? prob[:minloss] : -Inf
        if isfile(SOL_PATH)
            sol = open(SOL_PATH, "r") do f
                JSON.parse(f)
            end
        else
            sol = nothing
        end
        models = [
            # Symbol(Symbol("o$i") =>
            begin
                @unpack init, bbox = d
                L = bbox[2] - bbox[1]
                szd = Tuple(round.(Int, L / dx)) # design region size
                symmetries = [length(string(s)) == 1 ? Int(s) + 1 : s for s = d.symmetries]
                o = round((bbox[1] - origin) / dx) + 1
                # if !isa(init, AbstractArray)
                # if init == ""
                # if isnothing(init)
                #     init = eps_2D[o[1]:o[1]+szd[1]-1, o[2]:o[2]+szd[2]-1]
                #     # init = init .> (maximum(init) + minimum(init)) / 2
                #     init = init .≈ design_config.fill.epsilon |> F
                # elseif init == "random"
                #     init = nothing
                # end
                lmin = d.lmin / dx
                b = Blob(szd;
                    init=1,
                    lmin, rmin=lmin / 2,
                    contrast, T=F,
                    symmetries, verbose)
                if !isnothing(sol) && !restart
                    b.a .= sol.params[i] |> typeof(b.a)
                end
                b
            end for (i, d) = enumerate(designs)
        ]

        # model0 = deepcopy(model)
        # using CairoMakie
        # heatmap(model[1]())
    end

    boundaries = [] # unspecified boundaries default to PML

    device = 0
    device, dx, λc =
        F.((device, dx, λc))
    # guess = F.(guess)
    for ms = mode_solutions
        if !haskey(ms, :calibrated_modes)
            ms[:calibrated_modes] = []
        end
        for (mode, mode1) in zip(ms.modes, ms.modes1)

            sz = ms.size
            eps = ms.eps
            #  ϵmode1=   eps1 = ms.eps1

            mode = (; [k => complex.(stack.(v)...) |> F for (k, v) in mode |> pairs]...)
            mode1 = (; [k => complex.(v...) |> F for (k, v) in mode1 |> pairs]...)
            # ϵmode = eps |> stack |> transpose .|> F

            # sz = round(sz / dx) |> Tuple
            # # sz = size() - 1
            # mode = (; Pair.(keys(mode), [resize(v, size(v) - 1) for v in values(mode)])...)
            # ϵ = resize(ϵ, sz)
            if d == 2
                mode = mode1
                # mode = collapse_mode(mode,)
                # i = round(Int, size(ϵmode, 1) / 2)
                # ϵcore_ = ϵmode[i]
                # ϵmode = maximum(ϵmode, dims=2) |> vec
                # ϵmode = min.(ϵmode, ϵcore_)
            end
            # mode = normalize_mode(mode, dx / λc)
            mode = keepxy(mode)
            # global mode0 = deepcopy(mode)

            if calibrate && d == 2
                # @unpack mode, power = calibrate_mode(mode, ϵmode, dx / ms.wavelength)
                # global mode /= sqrt(power)
                # @unpack power, = calibrate_mode(mode, ϵmode, dx / λc;)
                # mode /= sqrt(power)

                # @unpack power, sol = calibrate_mode(mode, ϵmode, dx / λc; verbose=true)
                # MyMakie.save(joinpath(path, "calibration.png"), quickie(sol),)
                # @show power
                # global mode2 = deepcopy(mode)
            end

            # error()
            push!(ms[:calibrated_modes], mode)
        end
    end

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
                        λ = F(parse(F, string(λ)))
                        i = findfirst(mode_solutions) do v
                            isapprox(λ, v.wavelength) && string(port) in string.(v.ports)
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
                    λ = parse(F, string(λ))
                    λ / λc => [
                        begin
                            i = findfirst(mode_solutions) do v
                                isapprox(λ, v.wavelength) && port in v.ports && mn < length(v.modes)
                            end
                            ms = mode_solutions[i]
                            if isnothing(zcenter)
                                zcenter = ms.zcenter
                            end
                            mode = ms.calibrated_modes[mn+1]
                        end for mn = mn
                    ]
                end
                for (λ, mn) in pairs(m.wavelength_mode_numbers)
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
                ϵ = if run.d == 2
                    eps_2D
                    # heatmap(eps_2D) |> display
                else
                    eps_3D
                end
                sz = size(ϵ)
                ϵmax = maximum(ϵ)
                μ = 1
                σ = PML(1).σ
                δ = sqrt(ϵmax / μ) / σ
                setup(boundaries, sources, monitors, dx / λc, sz; F, Courant, ϵ,
                    pml_depths=[δ, δ, 0.3δ,],
                    pml_ramp_fracs=[0.2, 0.2, 1],
                    verbose,)
            end for (i, (run, sources, monitors)) in enumerate(zip(runs, runs_sources, runs_monitors))
        ]

    if !isempty(gpu_backend)
        println("using $gpu_backend backend.")
        Flux.gpu_backend!(gpu_backend)
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
    g0 = run_probs[1].geometry

    virgin = true
    # error()
    # Base.iterate(s::Symbol) = println(s)

    t0 = time()
    if study == "sparams"
        # @info "Computing s-parameters..."
        println("Computing s-parameters...")
        # sparams = write_sparams(img="", autodiff=false, verbose=true)
        sparams = write_sparams(runs, run_probs, g0, path, origin, dx;
            F, img="", verbose=true)
        sol = sparam_family(sparams)
    elseif study == "inverse_design"
        if length(origin) == 3
            if magic != "summersale"
                error("3D inverse design feature must be requested from Luminescent AI info@luminescentai.com")
            end
        end
        autodiff = true
        # compression = true
        global sparams = sparams0 = 0
        opt = Adam(eta)
        opt_state = Flux.setup(opt, models)
        println("starting optimization... first iter will be slow due to adjoint compilation.")
        stop = false
        global img = nothing
        best = best0 = 0
        S, ls = write_sparams(runs, run_probs, g0, path, origin, dx,
            designs, design_config, models;
            F, img, autodiff, with=true)
        # heatmap(_as[3])
        # global ass = gradient(models) do models
        #     write_sparams(runs, run_probs, g0, path, origin, dx,
        #         designs, design_config, models;
        #         F, img, autodiff, compression)
        # end
        # error()

        # prob = run_probs[1]
        # global aaaaa = gradient(g0) do geometry
        #     # solve(prob, geometry; autodiff, compression, verbose).forward_mode_powers[1][1][1]
        #     solve(prob, geometry; autodiff, compression, verbose)
        # end
        # error()

        for i = 1:iters
            # for i = 1:20
            # global virgin, stop, best, best0, sparams0
            global img = if virgin
                virgin = false
                "before.png"
            elseif i == iters || stop
                "after.png"
            end
            if isnothing(preset)
                @time l, (dldm,) = Flux.withgradient(models) do models
                    # sols = get_sols(runs, run_probs, g0, path, origin, dx,
                    S = write_sparams(runs, run_probs, g0, path, origin, dx,
                        designs, design_config, models;
                        F, img, autodiff, compression, ls)
                    l = 0
                    for k = keys(targets)
                        print("losses ")
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
                                            normalize(s1 / s2)
                                        end
                                    for ps = keys(targets[k][λ])])
                                for λ = keys(targets[k])])
                            y = fmap(cis, targets[k])
                            err = (x, y) -> angle(x / y) / F(π)
                        else
                            if :sparams == k
                                # S = get_sparams(sols)
                                ŷ = S
                            elseif :tparams == k
                                ŷ = fmap(abs2, S)
                            end
                            ŷ = [[ŷ(λ)(k) for k = keys(y[λ])] for λ = keys(y)]
                        end
                        _l = mean(abs, err.(flatten(ŷ), flatten(y)),) * weights(k)
                        print("$(k): $_l ")
                        l += _l
                    end
                    println("\n($i) weighted total loss $l")
                    l
                end
            elseif "phase_shifter" == preset.name
                @time l, (dldm,) = Flux.withgradient(models) do m
                    S = write_sparams(runs, run_probs, g0, path, origin, dx,
                        designs, design_config, m;
                        F, img, autodiff, compression, ls)#(1)(1)
                    k = keys(S) |> first
                    s = S[k][Symbol("o2@0,o1@0")]

                    S_ = write_sparams(runs, run_probs, g0, path, origin, dx,
                        designs, design_config, m;
                        F, img, autodiff, compression, perturb=:ϵ, ls)#(1)(1)
                    s_ = S_[k][Symbol("o2@0,o1@0")]

                    T = abs2(s)
                    dϕ = angle(s_ / s)
                    println("T: $T, dϕ: $dϕ")
                    global sparams = S
                    (T * dϕ)
                end
            end

            if i == 1
                best0 = best = l
            end

            Flux.update!(opt_state, models, dldm)# |> gpu)
            if l < best
                best = l
            end
            if l < minloss
                println("Loss below threshold, stopping optimization.")
                stop = true
                img = "after.png"
            end
            if i % 15 == 0
                if best - best0 > -0.01
                    println("Loss stagnating, stopping optimization.")
                    stop = true
                else
                    best0 = best
                end
            end
            println("")
        end
        println("Done in $(time() - t0) .")
        for (i, (m, d)) = enumerate(zip(models, designs))
            # Images.save(joinpath(path, "design$i.png"), Gray.(m() .< 0.5))
        end
        sol = (;
            sparam_family(sparams)...,
            optimized_designs=[m() .> 0.5 for m in models],
            params=getfield.(models, :a),
            designs,
            design_config,
        )
    end
    sol = (; sol..., path, dx, study,) |> cpu
    open(SOL_PATH, "w") do f
        write(f, json(sol))
    end
    sol
end