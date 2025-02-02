
function genrun(path, array=Array; kw...)
    Random.seed!(1)
    ENV["autodiff"] = "0"
    println("setting up simulation...")
    PROB = joinpath(path, "problem.json")
    SOL = joinpath(path, "solution.json")
    GEOMETRY = joinpath(path, "geometry")

    io = open(PROB)
    s = read(io, String)
    prob = JSON.parse(s; dicttype=OrderedDict)

    for (k, v) = kw
        prob[string(k)] = v
    end
    @unpack dtype, center_wavelength, dl, dx, xs, ys, zs, study, layer_stack, sources, monitors, materials, L, Ttrans, Tss, wavelengths = prob
    if study == "inverse_design"
        @unpack lsolid, lvoid, designs, targets, weights, eta, iters, restart, save_memory, design_config, stoploss = prob
    end

    ratio = int(dx / dl)
    F = Float32
    dtype = lowercase(dtype)
    if contains(dtype, "16") && contains(dtype, "bf")
        # F = BFloat16
        # println("BFloat16 selected. make sure your GPU supports it. otherwise will be emulated and very slow.")
    elseif contains(dtype, "16")
        F = Float16
        println("Float16 selected. make sure your cpu or GPU supports it. otherwise will be emulated and very slow.")
    end
    println("using $F")

    @show λ = center_wavelength
    @show λs = wavelengths / λ
    ticks = [
        begin
            round((v - v[1]) / dl)
        end for v = [xs, ys, zs]
    ]
    spacings = diff.(ticks)
    x = spacings[1][1]
    spacings = [x, x, spacings[3]]
    popfirst!.(ticks)
    deltas = spacings * dl
    mode_deltas = [dx, dx, dx]

    # global sz = round.(L / dl)

    GC.gc(true)
    if AUTODIFF()
        Nf = F
    else
        # Nf = [Bool, UInt8][findfirst(>=(length(enum), 2 .^ [1, 8]))]
        Nf = Float16
    end
    enum = [materials(v.material).epsilon |> F for v = values(layer_stack)] |> unique |> sort
    # ϵmin = minimum(enum) |> Nf

    _sz = Any[0]
    global geometry = DefaultDict(passkey=true) do k
        if k == :ϵ
            ones(Nf, _sz[1])
        elseif k == :σ
            zeros(Nf, _sz[1])
        end
    end


    # layer_stack = sort(collect(pairs(layer_stack)), by=kv -> -kv[2].mesh_order) |> OrderedDict
    # for (k, v) = pairs(layer_stack)
    ϵmin = 10000000
    for material = readdir(GEOMETRY)
        if string(material) ∉ ["sources", "monitors"]

            am = false
            for fn = readdir(joinpath(GEOMETRY, material))
                if endswith(lowercase(fn), "npy")
                    a = npzread(joinpath(GEOMETRY, material, "$fn"))
                    # a = permutedims(a, [2, 1, 3])
                    am = a .|| am
                end
            end
            _sz[1] = size(am)
            # println(k)
            # Figure()
            # display(volume(a))

            m = materials(material)
            if m(:epsilon) < ϵmin
                ϵmin = m(:epsilon)
            end
            for k = (:epsilon, :sigma)
                v = m(k, nothing)
                if !isnothing(v)
                    v = Nf(v)
                    @show k
                    k = togreek(k)
                    @show k
                    a = geometry[k]
                    a .*= .!(am)
                    a .+= am .* v
                end
            end
        end
    end

    geometry[:ϵ] = max.(ϵmin, geometry[:ϵ])
    GC.gc(true)
    geometry.ϵ |> volume |> display
    # error()


    global sources = map(enumerate(sources)) do (i, s)
        λsmode = (λs, s.mode)
        mask = npzread(joinpath(GEOMETRY, "sources", "$i.npy"))
        mask = permutedims(mask, [2, 1, 3])
        @assert !all(iszero, mask)
        mask = downsample(mask, ratio)
        Source(mask, F.(stack(s.frame)); λsmode)
    end
    global monitors = map(enumerate(monitors)) do (i, s)
        λsmode = (λs, s.mode)
        mask = npzread(joinpath(GEOMETRY, "monitors", "$i.npy"))
        mask = permutedims(mask, [2, 1, 3])
        @assert !all(iszero, mask)
        mask = downsample(mask, ratio)
        Monitor(mask, F.(stack(s.frame)); λsmode)
    end
    # volume(monitors[1].mask) |> display
    # error()

    boundaries = []
    N = 3
    @show Ttrans, Tss
    global prob = setup(dl / λ, boundaries, sources, monitors, deltas[1:N] / λ, mode_deltas[1:N-1] / λ;
        geometry..., array, F, deltas3=deltas / λ, Ttrans, Tss,
        # lpml=[0.2, 0.2, 0.2],
        lpml=ones(3),)#
    # σpml=4,)#
    # pml_depths=[0.2, 0.2, 0.2])

    # v = prob.monitor_instances
    # v[1].frame = nothing
    # v[1].dimsperm = [1, 2, 3]
    # v[2].frame = nothing
    # v[2].dimsperm = [-1, 2, -3]

    global sol = solve(prob; path)
    global S = getsparams.((sol,), eachindex(λs))

    # solution = (;
    #     sparam_family(S)...,)
    solution = sparam_family(S)
    json(solution, 4) |> println
    # volume(prob._geometry.σ |> cpu)
    # volume(geometry.σ |> cpu)
    open(joinpath(path, "solution.json"), "w") do f
        write(f, json(cpu(solution)))
    end
    u = cpu(sol.u)
    d = merge(u.E, u.H)
    d = kmap(string, identity, d) |> pairs |> Dict
    npzwrite(joinpath(path, "fields.npz"), d)

    plotslices(prob._geometry.ϵ, joinpath(path, "epsilon.png"))
    plotslices(u.Ey, joinpath(path, "Ey.png"))
    S, sol, prob
end