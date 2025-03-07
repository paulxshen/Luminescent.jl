
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
    @unpack dorun, dtype, center_wavelength, dl, dx, xs, ys, zs, study, layer_stack, sources, monitors, materials, L, Ttrans, Tss, Tssmin, wavelengths, frequencies, df = prob
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
    plane_deltas = [dx, dx, dx]

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
    # ϵmin = 10000000
    for material = readdir(GEOMETRY)
        if string(material) ∉ ["sources", "monitors"]
            @show material
            am = false
            for fn = readdir(joinpath(GEOMETRY, material))
                if endswith(lowercase(fn), "npy")
                    a = npzread(joinpath(GEOMETRY, material, "$fn"))
                    # a = permutedims(a, [2, 1, 3])
                    am = a .|| am
                end
            end
            _sz[1] = size(am)
            # display(volume(am))
            # error()

            m = materials(material)
            # if m(:epsilon) < ϵmin
            # ϵmin = m(:epsilon)
            # end
            for k = (:epsilon, :sigma)
                v = m(k, nothing)
                if !isnothing(v)
                    @show k
                    v = Nf(v)
                    k = togreek(k)
                    a = geometry[k]
                    a .*= .!(am)
                    a .+= am .* v
                end
            end
        end
    end

    ϵmin = 1
    geometry[:ϵ] = max.(ϵmin, geometry[:ϵ])
    GC.gc(true)
    # geometry.ϵ |> volume |> display
    # error()


    global sources = map(enumerate(sources)) do (i, s)
        λsmode = (λs, s.mode)
        mask = npzread(joinpath(GEOMETRY, "sources", "$i.npy"))
        @assert !all(iszero, mask)
        mask = downsample(mask, ratio)
        Source(mask, F.(stack(s.frame)); λsmode)
    end
    global monitors = map(enumerate(monitors)) do (i, s)
        λsmode = (λs, s.mode)
        mask = npzread(joinpath(GEOMETRY, "monitors", "$i.npy"))
        @assert !all(iszero, mask)
        mask = downsample(mask, ratio)
        Monitor(mask, F.(stack(s.frame)); λsmode)
    end
    # volume(sources[1].mask) |> display
    # error()

    boundaries = []
    N = 3
    global prob = setup(dl / λ, boundaries, sources, monitors, deltas[1:N] / λ, plane_deltas[1:N-1] / λ;
        geometry..., array, F, deltas3=deltas / λ, Ttrans, Tss, Tssmin,
        # pmlfracs=[0.2, 0.2, 0.2],
        # pmlfracs=1 / λs[1] * ones(3),
        # σpml=4,#
        # pml_depths=[0.2, 0.2, 0.2])
    )
    global g = prob._geometry.ϵ |> cpu
    g = min.(g, 70)
    # volume(g) |> display
    # error()
    plotslices(g; saturation=10, path=joinpath(path, "epsilon.png"))
    !dorun && return prob

    # v = prob.monitor_instances
    # v[1].frame = nothing
    # v[1].dimsperm = [1, 2, 3]
    # v[2].frame = nothing
    # v[2].dimsperm = [-1, 2, -3]

    global sol = solve(prob; path)

    fs = range(frequencies[1], 1.0001frequencies[end], step=df)
    global S = getsparams.((sol,), reverse(eachindex(λs)))

    sarray = vec2smatrix.(S)
    Is = CartesianIndices(sarray[1])
    frequencies = range(frequencies[1], frequencies[end], length=length(frequencies))
    sarray = map(Is) do I
        v = getindex.(sarray, I)
        itp = cubic_spline_interpolation(frequencies, v)
        itp.(fs)
    end
    sarray = stack(map(eachindex(sarray[1])) do i
        getindex.(sarray, i)
    end)
    sarray = permutedims(sarray, [3, 1, 2])

    sarray = ComplexF32.(sarray)
    npzwrite(joinpath(path, "S.npy"), sarray)

    # solution = (;
    #     sparam_family(S)...,)
    solution = sparam_family(S)
    json(solution, 4) |> println
    # volume(prob._geometry.σ |> cpu)
    # volume(geometry.σ |> cpu)
    solution[:fs] = fs
    open(joinpath(path, "solution.json"), "w") do f
        write(f, json(cpu(solution)))
    end

    u = cpu(sol.u)
    d = merge(u.E, u.H)
    global d = kmap(string, identity, d) |> pairs |> Dict
    # npzwrite(joinpath(path, "fields.npz"), d)


    plotslices(d.Ey |> cpu; saturation=1e3, path=joinpath(path, "Ey.png"))
    # for saturation = [1, 5]
    # end

    # plotslices(d.Ey |> cpu; saturation=1)
    # volume(sol.u.E.Ey |> cpu) |> display

    # c = sol.u.E.Ey
    # sz = c |> size
    # a = zeros(prob.grid.sz)
    # @assert sz == size(a)
    # setindexf!(a, 1, prob.monitor_instances[1].inds(1)...)
    # setindexf!(a, 2, prob.monitor_instances[2].inds(1)...)
    # b = abs.(prob.source_instances[1].sigmodes(1)[2].Jy) |> cpu
    # b /= maximum(b)
    # a .+= 3b
    # volume(a) |> display

    S, sol, prob
end