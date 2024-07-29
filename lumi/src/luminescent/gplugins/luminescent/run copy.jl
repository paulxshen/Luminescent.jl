#=
We do inverse design of a compact photonic waveguide bend to demonstrate workflow of FDTD adjoint optimization. First, we seed the design using 2d TE adjoint simulations which serve as fast approximations. Optionlly, we finetune the resulting design in full blown 3d adjoint simulations.

=#

using UnPack, LinearAlgebra, Random, StatsBase, Dates, Colors, DataStructures, JSON, Images
using Zygote, Flux, GLMakie, Jello, ImageTransformations
using Flux: mae, Adam
using Zygote: withgradient, Buffer
using BSON: @save, @load
using AbbreviatedStackTraces
using Rsvg, Cairo, FileIO
# using Jello, Luminescent, LuminescentVisualization
# using Luminescent:keys, values
Random.seed!(1)

# if running directly without module # hide
include("$(pwd())/src/main.jl") # hide
include("$(pwd())/../LuminescentVisualization.jl/src/main.jl") # hide
# using Porcupine: keys, values

if isempty(ARGS)
    path = joinpath(pwd(), "lumi_runs")
    path = filter(isdir, readdir(path, join=true)) |> sort |> last
    PROB_PATH = joinpath(path, "prob.bson")
else
    PROB_PATH = ARGS[1]
    path = dirname(PROB_PATH)

end
#=
We skip 3d finetuning as it's 20x more compute and memory intensive than 2d adjoints. If wishing to do 3d finetuning, set `iterations3d`. In any case, 3d forward simulations (without adjoint) only take a few seconds.
=#
# name = "ubend"
calibrate = false
F = Float32
model_name = nothing # if load saved model
tol = 1e-3
dosave = false
verbose = false
#=
We load design layout which includes a 2d device of device waveguide geometry as well as variables with locations of ports, sources, design regions and material properties.
=#
@load PROB_PATH name runs ports λc dx components study mode_solutions eps_2D eps_3D mode_height zmin gpu
eps_2D = F(stack(stack.(eps_2D)))'
eps_3D = F(permutedims(stack(stack.(eps_3D)), (3, 2, 1)))
# heatmap(eps_2D) |> display
# GLMakie.volume(eps_3D) |> display
# error()
#  ϵbase ϵclad ϵcore hbase hwg hclad
ports, = (ports,) .|> SortedDict
polarization = :TE

_dx = dx / λc
m = round(SOURCE_MARGIN / dx)
source_margin = m * dx
n = round(PORT_SOURCE_OFFSET / dx)
port_source_offset = n * dx
# device = pad(device, :replicate, n)
eps_2D = pad(eps_2D, :replicate, m + n)
eps_3D = pad(eps_3D, :replicate, (m + n, m + n, 0))
origin = components.device.bbox[1] - (m + n) * dx

if study == "inverse_design"
    @load PROB_PATH designs design_region_layer targets maxiters design_config
    #=
    We initialize a Jello.jl Blob object which will generate geometry of design region. Its parameters will get optimized during adjoint optimization. We initialize it with a straight slab connecting input to output port.
    =#

    model = NamedTuple([
        Symbol("m$i") => begin
            bbox = d.bbox
            L = bbox[2] - bbox[1]
            szd = Tuple(round.(Int, L / dx)) # design region size
            s = d.symmetries
            o = round((bbox[1] - origin) / dx) + 1
            init = 1
            if !isa(init, AbstractArray)
                init = eps_2D[o[1]:o[1]+szd[1]-1, o[2]:o[2]+szd[2]-1]
                # init = init .> (maximum(init) + minimum(init)) / 2
                global init = init .≈ design_config.fill.ϵ |> F
            end
            Blob(szd;
                init,
                lmin=d.lmin / dx,
                contrast=10,
                rmin=nothing,
                symmetry_dims=isempty(s) ? s : s + 1)
        end for (i, d) = enumerate(designs)
    ])
    model = model.m1

    targets = SortedDict([F(λ) => SortedDict([(k) => v for (k, v) = pairs(d)]) for (λ, d) = pairs(targets)])
    # model0 = deepcopy(model)
    # heatmap(model())
end
# Flux.trainable(v::AbstractArray{Blob}) = NamedTuple([Symbol("a$i") => Flux.trainable(model) for (i, model) = enumerate(v)])
# Flux.trainable(model)
# error()
#=
We set boundary conditions, sources , and monitor. The modal source profile is obtained from external mode solver , in our case VectorModesolver.jl . Please refer to guide section of docs website for details . To get an approximate  line source for use in 2d from the cross section profile , we sum and collapse it along its height axis
=#

boundaries = [] # unspecified boundaries default to PML

# monitors = [monitors[1]]

# device, guess,  Δ, ϵbase, ϵcore, ϵclad, dx, λc =
#     F.((device, guess,  Δ, ϵbase, ϵcore, ϵclad, dx, λc))

device = 0
device, dx, λc =
    F.((device, dx, λc))
# guess = F.(guess)

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
                    if !haskey(ms, :calibrated_modes)
                        ms[:calibrated_modes] = Any[nothing for i = eachindex(ms.modes)]
                    end
                    mode = ms[:calibrated_modes][mn+1]
                    if isnothing(mode)
                        mode = ms.modes[mn+1]

                        sz = ms.size
                        eps = ms.eps

                        mode = (; [k => transpose(complex.(stack.(v)...)) |> F for (k, v) in mode |> pairs]...)
                        ϵmode = eps |> stack |> transpose .|> F

                        # sz = round(sz / dx) |> Tuple
                        # # sz = size() - 1
                        # mode = (; Pair.(keys(mode), [resize(v, size(v) - 1) for v in values(mode)])...)
                        # ϵ = resize(ϵ, sz)
                        if d == 2
                            mode, ϵmode = collapse_mode(mode, ϵmode)
                            i = round(Int, size(ϵmode, 1) / 2)
                            ϵcore_ = ϵmode[i]
                            ϵmode = maximum(ϵmode, dims=2) |> vec
                            ϵmode = min.(ϵmode, ϵcore_)
                        end
                        # mode = normalize_mode(mode, dx / λc)
                        mode = keepxy(mode)
                        global mode0 = deepcopy(mode)

                        if calibrate
                            @unpack mode, power = calibrate_mode(mode, ϵmode, dx / λc)
                            # global mode /= sqrt(power)
                            # plot(abs.(mode.Ex)) |> display
                            # plot(abs.(mode0.Ex)) |> display
                            @unpack power, = calibrate_mode(mode, ϵmode, dx / λc;)
                            mode /= sqrt(power)

                            @unpack power, sol = calibrate_mode(mode, ϵmode, dx / λc; verbose=true)
                            GLMakie.save(joinpath(path, "calibration.png"), quickie(sol),)
                            @show power
                            # global mode2 = deepcopy(mode)
                        end

                        # error()
                        ms[:calibrated_modes][mn+1] = mode
                    end
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
        wavelength_modes = OrderedDict([
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
        prob = setup(boundaries, sources, monitors, dx / λc, sz; F, ϵ)
        global bb = prob
    end for (i, (run, sources, monitors)) in enumerate(zip(runs, runs_sources, runs_monitors))
]

if gpu
    CUDA.allowscalar(false)
    @assert CUDA.functional()
    run_probs = gpu.(run_probs)
end

function write_sparams(model=nothing)
    sol = [
        begin

            if !isnothing(model)
                # prob[:geometry] = make_geometry(model)
                # geometry = make_geometry(model)
                # sol = solve(prob, geometry; verbose,)
                sol = solve(prob, model; verbose)
            else
                sol = solve(prob; verbose)
            end

            ignore() do
                try
                    GLMakie.save(joinpath(path, "run$i.png"), quickie(sol),)
                catch e
                    @show e
                end
            end
            sol
        end for (i, prob) in enumerate(run_probs)
    ]

    coeffs = OrderedDict()
    for (v, run) in zip(sol, runs)
        sources = run.sources |> values
        monitors = run.monitors |> values
        source_port = first(values(sources)).port
        source_mn = values(first(values(sources)).wavelength_mode_numbers)[1][1]
        for (monitor, v) = zip(monitors, v.mode_coeffs)
            for (λ, v) = zip(monitor.wavelength_mode_numbers |> keys, v)
                for (monitor_mn, v) = zip(monitor.wavelength_mode_numbers[λ], v)
                    monitor_port = monitor.port
                    if !haskey(coeffs, λ)
                        coeffs[λ] = OrderedDict()
                    end
                    coeffs[λ]["o$monitor_port@$monitor_mn,o$source_port@$source_mn"] = v
                end
            end
        end
    end
    sparams = dict([λ => dict([k => begin
        s = 0
        ignore() do
            s = split(k, ",")[2]
        end
        # k = string(k)
        # i = findfirst(',', string(k))
        # s = string(k)[i+1:length(k)]
        # s = "o1@0"
        coeffs[λ][k][1] / coeffs[λ]["$s,$s"][2]
    end for (k) = keys(coeffs[λ])]) for λ = keys(coeffs)])

    # sparams = OrderedDict([λ => OrderedDict([k => begin
    #     s = split(k, ",")[2]
    #     v[1] / d["$s,$s"][2]
    # end for (k, v) = pairs(d)]) for (λ, d) = pairs(coeffs)])
    # if source_mn == monitor_mn == 0
    #     coeffs[λ]["$monitor_port,$source_port"] = v
    # end
    sparams
end

g0 = run_probs[1].geometry |> deepcopy
# props = keys(g0)
function make_geometry(model)
    # dict([k => begin
    #     a = g0[k]
    #     # a = ones(F, size(g0[k]))
    #     if haskey(design_config.fill, k)
    #         for (m, d) in zip(model, designs)
    #             mask = m()

    #             println(k)
    #             # if k in props

    #             f = design_config.fill[k]
    #             v = design_config.void[k]
    #             # for (f, v) = zip(design_config.fill, design_config.void)
    #             p = (v .* (1 - mask) + f .* mask)
    #             a = place(a, ((d.bbox[1] - origin) / dx) + 1, p; replace=true) |> F
    #             # a = g0[k]
    #             # b = Zygote.Buffer(a)
    #             # copyto!(b, a)
    #             # place!(b, round((d.bbox[1] - origin) / dx) + 1, p |> F)
    #             # g[k] = copy(b)
    #         end
    #     end
    #     a
    # end for k = keys(g0,)])
    mask = model[1]()
    f = design_config.fill.ϵ
    v = design_config.void.ϵ
    # for (f, v) = zip(design_config.fill, design_config.void)
    p = (v .* (1 - mask) + f .* mask)
    a = place(g0.ϵ, ((designs[1].bbox[1] - origin) / dx) + 1, p; replace=true) |> F
    (; ϵ=a, μ=ones(F, size(a)), σ=ones(F, size(a)), σm=ones(F, size(a)))
end
# g = make_geometry(model)
# f, a, pl = heatmap(g[:ϵ])
# f, a, pl = heatmap(p)
# Colorbar(f[1, 2], pl)
# f

function loss(sparams,)# targets)
    # prob[:geometry] = make_geometry(device, ϵ1, ϵ2, model, model_origin)
    mean([mae(abs.(getindex.((yhat,), string.(keys(y)))), values(y).mag) for (y, yhat) = zip(values(targets), values(sparams))])
    # mean([mae(values(yhat),values(y)) for (yhat,y)=zip(targets,[solve(prob;geometry) for prob=run_probs])])
end
function sparam_family(sparams)

    tparams = apply(abs2, sparams)
    sparam_phasors = apply(sparams) do z
        (abs(z), angle(z))
    end
    println(tparams)
    # sparams = apply(reim, sparams)
    sol = (; sparams, tparams, sparam_phasors)
end
if study == "sparams"
    sparams = write_sparams()
    # tparams = dict([k => abs(sparams[k])^2 for k = keys(sparams)])
    # sol = [(; wavelength=λ, sparams=dict([k => reim(v) for (k, v) = pairs(d)])) for (λ, d) = pairs(sparams)]
    sol = sparam_family(sparams)
elseif study == "inverse_design"
    model = model
    # @show loss(model, targets)
    opt = Adam(1)
    opt_state = Flux.setup(opt, model)
    # for i = 1:maxiters
    for i = 1:5
        println("$i")
        @time global l, (dldm,) = withgradient(model) do m
            sparams = write_sparams([m])
            loss(sparams)
        end
        l < 0 && break
        Flux.update!(opt_state, model, dldm)
        println(" $l\n")
    end
    for (i, (m, d)) = enumerate(zip(model, designs))
        Images.save(joinpath(path, "design$i.png"), Gray.(m() .< 0.5))
    end
    model = model
    sol = (;
        sparam_family(sparams)...,
        designs,

        # model,
        # designs=[m() for m in model],
    )
end
# @save "$path/sol.json" sol
open("$(path)/sol.json", "w") do f
    write(f, json(sol))
end