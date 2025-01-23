function genrun(path)
    PROB = joinpath(path, "problem.json")
    prob = JSON.parse(read(open(PROB), String); dicttype=OrderedDict)
    gpu_backend = prob["gpu_backend"]
    array = if isnothing(gpu_backend)
        println("using CPU")
        Array
    else
        println("using $gpu_backend")
        cu
    end
    picrun(path, array)
end

function genrun(path, array; kw...)
    Random.seed!(1)
    ENV["autodiff"] = "0"
    println("setting up simulation...")
    PROB = joinpath(path, "problem.json")
    SOL = joinpath(path, "solution.json")
    TEMP = joinpath(path, "TEMP")

    io = open(PROB)
    s = read(io, String)
    global prob = JSON.parse(s; dicttype=OrderedDict)


    @unpack dtype, center_wavelength, dl, xs, ys, zs, study, layer_stack, sources, monitors, materials, L, Ttrans, Tss, wavelengths = prob
    if study == "inverse_design"
        @unpack lsolid, lvoid, designs, targets, weights, eta, iters, restart, save_memory, design_config, stoploss = prob
    end

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
    λs = wavelengths / λ
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


    # global sz = round.(L / dl)

    GC.gc(true)
    if AUTODIFF()
        Nf = F
    else
        # Nf = [Bool, UInt8][findfirst(>=(length(enum), 2 .^ [1, 8]))]
        Nf = Float16
    end
    enum = [materials(v.material).epsilon |> F for v = values(layer_stack)] |> unique |> sort
    ϵmin = minimum(enum) |> Nf
    ϵ3 = nothing

    # layer_stack = sort(collect(pairs(layer_stack)), by=kv -> -kv[2].mesh_order) |> OrderedDict
    for (k, v) = pairs(layer_stack)
        @unpack material = v
        ϵ = materials(material).epsilon |> Nf
        a = npzread(joinpath(TEMP, "$k.npy"))
        a = permutedims(a, [2, 1, 3])

        if isnothing(ϵ3)
            ϵ3 = zeros(Nf, size(a))
        end
        ϵ3 .*= 1 - a
        ϵ3 .+= a .* ϵ
    end
    ϵ3 = max.(ϵmin, ϵ3)
    @assert eltype(ϵ3) == Nf
    GC.gc(true)

    sources = map(enumerate(sources)) do (i, s)
        λsmode = (λs, s.mode)
        mask = npzread(joinpath(TEMP, "source$i.npy"))
        Source(mask; λsmode)
    end
    monitors = map(enumerate(monitors)) do (i, s)
        λsmode = (λs, s.mode)
        mask = npzread(joinpath(TEMP, "monitor$i.npy"))
        Monitor(mask; λsmode)
    end

    boundaries = []
    ϵ = ϵ3
    N = 3
    prob = setup(dl / λ, boundaries, sources, monitors, deltas[1:N] / λ, mode_deltas[1:N-1] / λ, ; array,
        F, ϵ, deltas3=deltas / λ, λ, TEMP, Ttrans, Tss)


    sol = solve(prob; path)
end