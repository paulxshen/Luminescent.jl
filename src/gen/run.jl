
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
    global ϵ3 = nothing

    # layer_stack = sort(collect(pairs(layer_stack)), by=kv -> -kv[2].mesh_order) |> OrderedDict
    # for (k, v) = pairs(layer_stack)
    for material = readdir(GEOMETRY)
        if string(material) ∉ ["sources", "monitors"]
            ϵ = materials(material).epsilon |> Nf
            am = false
            for fn = readdir(joinpath(GEOMETRY, material))
                if endswith(lowercase(fn), "npy")
                    a = npzread(joinpath(GEOMETRY, material, "$fn"))
                    a = permutedims(a, [2, 1, 3])
                    am = a .|| am
                end
            end
            # println(k)
            # Figure()
            # display(volume(a))

            if isnothing(ϵ3)
                ϵ3 = ones(Nf, size(am))
            end
            ϵ3 .*= .!(am)
            ϵ3 .+= am .* ϵ
        end
    end

    # ϵ3 = max.(ϵmin, ϵ3)
    @assert eltype(ϵ3) == Nf
    GC.gc(true)

    # volume(ϵ3) |> display
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
    # error()

    boundaries = []
    ϵ = ϵ3
    N = 3
    @show Ttrans, Tss
    global prob = setup(dl / λ, boundaries, sources, monitors, deltas[1:N] / λ, mode_deltas[1:N-1] / λ;
        Courant=0.5, array, F, ϵ, deltas3=deltas / λ, Ttrans, Tss, lpml=[0.3, 0.3, 0.5],)# σpml=0)


    sol = solve(prob; path)
    prob, sol
end