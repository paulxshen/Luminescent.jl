
using Luminescent, Random

Random.seed!(1)
if isempty(ARGS)
    path = lastrun()
    println("path: ", path)
else
    path = ARGS[1]
end
# gfrun(path)

using BSON: @load, @save, load
using BSON, Statistics, DataStructures, Porcupine, ArrayPadding, UnPack, Jello, JSON, Flux, CUDA, Zygote
using Luminescent: gpu, cpu
println("setting up simulation...")
# do something based on ARGS?
PROB_PATH = joinpath(path, "prob.bson")

calibrate = true
model_name = nothing # if load saved model
tol = 1e-3
dosave = false
verbose = false
#=
We load design layout which includes a 2d device of device waveguide geometry as well as variables with locations of ports, sources, design regions and material properties.
=#
@load PROB_PATH name dtype margin source_margin Courant port_source_offset portsides runs ports dx components study mode_solutions eps_2D eps_3D mode_height zmin gpu_backend d
F = Float32
if contains(dtype, "16")
    # F = Float16
    println("Float16 not supported yet, will be in future release.")
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
source_margin = whole(source_margin, dx)

p = round((port_source_offset + source_margin) / dx)
np = round(margin / dx)
# p = m + n
# origin = components.device.bbox[1] - dx * p
# eps_3D = pad(eps_3D, :replicate, [p, p, 0])
# eps_2D = pad(eps_2D, :replicate, p)
model = nothing
lr = p * portsides + np * (1 - portsides)
eps_2D = pad(eps_2D, :replicate, lr...)
# nz=whole(ZMARGIN,_dx)
push!(lr[1], 0)
push!(lr[2], 0)
eps_3D = pad(eps_3D, :replicate, lr...)
# heatmap(eps_2D) |> display
origin = components.device.bbox[1] - dx * (p * portsides[1] + np * (1 - portsides[1]))

if study == "inverse_design"
    @load PROB_PATH designs design_layer targets target_type eta maxiters design_config
    prob = load(PROB_PATH)
    minloss = haskey(prob, :minloss) ? prob[:minloss] : -Inf
    #=
    We initialize a Jello.jl Blob object which will generate geometry of design region. Its parameters will get optimized during adjoint optimization. We initialize it with a straight slab connecting input to output port.
    =#

    model = [
        # Symbol("m$i") =>
        begin
            @unpack init, bbox = d
            L = bbox[2] - bbox[1]
            szd = Tuple(round.(Int, L / dx)) # design region size
            symmetry_dims = [length(string(s)) == 1 ? Int(s) + 1 : s for s = d.symmetries]
            o = round((bbox[1] - origin) / dx) + 1
            # if !isa(init, AbstractArray)
            # if init == ""
            if isnothing(init)
                init = eps_2D[o[1]:o[1]+szd[1]-1, o[2]:o[2]+szd[2]-1]
                # init = init .> (maximum(init) + minimum(init)) / 2
                init = init .≈ design_config.fill.ϵ |> F
            elseif init == "random"
                init = nothing
            end
            lmin = d.lmin / dx
            Blob(szd;
                init, lmin, rmin=lmin / 2,
                contrast=10,
                symmetry_dims, verbose)
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
            # mode, ϵmode = collapse_mode(mode, ϵmode)
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
                    L = [sig.width] / λc
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
        L = [m.width] / λc

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

run_probs = [
    begin
        ϵ = if run.d == 2
            eps_2D
            # heatmap(eps_2D) |> display
        else
            eps_3D
        end
        sz = size(ϵ)
        prob = setup(boundaries, sources, monitors, dx / λc, sz; F, Courant, ϵ, zpml=0.2, verbose)
    end for (i, (run, sources, monitors)) in enumerate(zip(runs, runs_sources, runs_monitors))
]
# error()
if !isempty(gpu_backend)
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
    model = model |> gpu
end
g0 = run_probs[1].geometry |> deepcopy
virgin = true

t0 = time()
if study == "sparams"
    # @info "Computing s-parameters..."
    println("Computing s-parameters...")
    # sparams = write_sparams(img="", autodiff=false, verbose=true)
    sparams = write_sparams(runs, run_probs, g0, path, origin, dx;
        F, img="", autodiff=true, verbose=true)
    sol = sparam_family(sparams)
elseif study == "inverse_design"
    sparams0 = 0
    # @show write_sparams(model)
    opt = Adam(eta)
    opt_state = Flux.setup(opt, model)
    # @info "starting optimization... first iter will be slow due to compilation."
    println("starting optimization... first iter will be slow due to adjoint compilation.")
    stop = false
    img = nothing
    best = Inf
    global virgin = true
    for i = 1:maxiters
        img = if virgin
            global virgin = false
            "before.png"
        elseif i == maxiters || stop
            "after.png"
        else
            img
        end
        @time l, (dldm,) = withgradient(model) do m
            sparams = write_sparams(runs, run_probs, g0, path, origin, dx,
                designs, design_config, m, ;
                F, img, autodiff=true)
            params = if target_type == "sparams"
                sparams
            else
                Porcupine.apply(abs2, sparams)
            end
            loss(params, targets)
            # abs(sparams[1].mode_coeffs[2][1][1][1])

        end
        if stop
            break
        end
        if i == 1
            best0 = best = l
            sparams0 = deepcopy(sparams)
        end
        # l < 0 && break
        Flux.update!(opt_state, model, dldm)
        if l < best
            best = l
        end
        if l < minloss
            # @info "Loss below threshold, stopping optimization."
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
        println("$i loss $l\n")
    end
    # @info "Done in $(time() - t0) ."
    println("Done in $(time() - t0) .")
    for (i, (m, d)) = enumerate(zip(model, designs))
        Images.save(joinpath(path, "design$i.png"), Gray.(m() .< 0.5))
    end
    sol = (;
        before=sparam_family(sparams0),
        after=sparam_family(sparams),
        design_config, path, dx,
        designs=[m() .> 0.5 for m in model],
    )
end
# @save "$path/sol.json" sol
open("$(path)/sol.json", "w") do f
    write(f, json(sol))
end
sol