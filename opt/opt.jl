#=
We do inverse design of a compact photonic waveguide bend to demonstrate workflow of FDTD adjoint optimization. First, we seed the design using 2d TE adjoint simulations which serve as fast approximations. Optionlly, we finetune the resulting design in full blown 3d adjoint simulations.

=#

using UnPack, LinearAlgebra, Random, StatsBase, Dates, Colors, DataStructures
using Zygote, Flux, CUDA, GLMakie, Jello, Images
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

#=
We skip 3d finetuning as it's 20x more compute and memory intensive than 2d adjoints. If wishing to do 3d finetuning, set `iterations3d`. In any case, 3d forward simulations (without adjoint) only take a few seconds.
=#
# name = "ubend"
iterations2d = 0
record2d = true
record3d = false
F = Float32
ongpu = false
model_name = nothing # if load saved model
tol = 1e-3
dosave = false

#=
We load design layout which includes a 2d device of device waveguide geometry as well as variables with locations of ports, signals, design regions and material properties.
=#
@load "$(@__DIR__)/prob.bson" name signals ports designs λ dx ϵbase ϵclad ϵcore hbase hwg hclad layers path_length targets margins lmin
signals, ports, targets = (signals, ports, targets) .|> SortedDict
polarization = :TE

# error()
for s = [:device, :init]
    # for name = ["device", "init"]
    r = Rsvg.handle_new_from_data(layers[s][:svg])
    # r = Rsvg.handle_new_from_file("$(@__DIR__)/$name.svg")
    d = Rsvg.handle_get_dimensions(r)
    cs = Cairo.CairoImageSurface(1 * (d.width), 1 * (d.height), Cairo.FORMAT_ARGB32)
    c = Cairo.CairoContext(cs)
    Rsvg.handle_render_cairo(c, r)
    Cairo.write_to_png(cs, "$(@__DIR__)/$s.png")
end

masks = Dict([k => imresize(F.(convert.(Gray, FileIO.load("$(@__DIR__)/$k.png") |> transpose),), round.(Int, (layers[k].bbox[2] - layers[k].bbox[1]) / dx) |> Tuple) .> tol for k in [:device, :init]])
@unpack device, init = masks
# heatmap(masks.device) |> display
lm, rm = round.(margins / dx)
device = pad(device, 0, lm, rm)
origin = layers.device.bbox[1] - lm * dx
#=
We initialize a Jello.jl Blob object which will generate geometry of design region. Its parameters will get optimized during adjoint optimization. We initialize it with a straight slab connecting input to output port.
=#

bbox = designs[1].bbox
L = bbox[2] - bbox[1]
szd = Tuple(round.(Int, L / dx)) # design region size
model = Blob(szd; init,
    lmin=lmin / dx,
    contrast=10,
    rmin=nothing,
    symmetry_dims=(designs[1].symmetry_dims.|>Int)[1])
model0 = deepcopy(model)
heatmap(model())

#=
We set key time intervals. The signal must first propagate to port 2 after which all port power fluxes will get monitored
=#

Δ = zeros(2)
# Δ[1] = 1
Δ[1] = 2 + path_length / λ * sqrt(ϵcore) # simulation duration in [periods] for signal to reach output ports
Δ[2] = 2 # duration to record power at output ports
T = cumsum(Δ)

#=
We set boundary conditions, sources , and monitor. The modal source profile is obtained from external mode solver , in our case VectorModesolver.jl . Please refer to guide section of docs website for details . To get an approximate  line source for use in 2d from the cross section profile , we sum and collapse it along its height axis
=#

boundaries = [] # unspecified boundaries default to PML
sig = (signals |> values |> first)
wm = sig.width
monitors = [
    # (center, lower bound, upper bound; normal)
    begin
        c = (p.center - origin) / λ
        n = p.normal
        if k in keys(signals)
            c -= margin_offset * n
            c -= source_monitor_offset * n
            n *= -1
        else
            c -= source_monitor_offset * n
        end
        L = (p.endpoints[2] - p.endpoints[1])
        L = L / norm(L) * wm / λ
        Monitor(c, L, n)
    end
    for (k, p) = ports |> pairs
]
# monitors = [monitors[1]]

# modal source
mode = (; [k => transpose(complex.(stack.(v)...)) for (k, v) in sig.mode |> pairs]...)
ϵ = sig.ϵ |> stack |> transpose .|> F

sz = round(sig.size / dx) |> Tuple
mode = (; Pair.(keys(mode), [resize(v, sz) for v in values(mode)])...)
ϵ = resize(ϵ, sz)

mode, ϵmode = collapse_mode(mode, ϵ)
mode = normalize_mode(mode, dx / λ)
i = round(Int, size(ϵmode, 1) / 2)
ϵcore_ = ϵmode[i]
ϵslice = maximum(ϵ, dims=2) |> vec
ϵslice = min.(ϵslice, ϵcore_)

mode, a, b = calibrate_mode(mode, ϵslice, dx / λ; name="1")
mode = normalize_mode(mode, dx / λ)
mode_, mp0, = calibrate_mode(mode, ϵslice, dx / λ; name="2")
@show mp0
# @unpack Ex, Ez = mode

# error()

# zaxis = sig.normal
@unpack center, endpoints, = sig
n = sig.normal
c = (sig.center - origin) / λ + margin_offset * sig.normal
L = (sig.endpoints[2] - sig.endpoints[1]) / λ
# sources = [ModalSource(t -> cispi(2t), c, L, n; Jx=mode.Ex, Jz=mode.Ez)]
sources = [ModalSource(t -> cispi(2t), c, L, n; Jx=mode.Ex,)]

ϵmin = ϵclad
sz = size(device)

device, init, T, Δ, ϵbase, ϵcore, ϵcore_, ϵclad, dx, mp0, =
    F.((device, init, T, Δ, ϵbase, ϵcore, ϵcore_, ϵclad, dx, mp0,))

configs = maxwell_setup(boundaries, sources, monitors, dx / λ, sz; F, ϵmin, polarization,)
@unpack dt, sz, geometry_padding, geometry_staggering, field_padding, source_instances, monitor_instances, u0, polarization = configs
nt = round(Int, 1 / dt)
A = area.(monitor_instances)
# nmp = nfp.Ex ⋅ nfp.Hy * A[1] / length(Ex) / 2

if ongpu
    using Flux
    # using CUDA
    # @assert CUDA.functional()
    u0, model, device, μ, σ, σm, field_padding, source_instances =
        gpu.((u0, model, device, μ, σ, σm, field_padding, source_instances))
    merge!(configs, (; u0, field_padding, source_instances))
end

#=
We define a geometry update function that'll be called each adjoint iteration. It calls geometry generator model to generate design region which gets placed onto mask of device features.
    =#

function make_geometry(model, configs, dϵ=0)#; make3d=false)
    @unpack sz, geometry_padding, geometry_staggering = configs
    μ = ones(F, sz)
    σ = zeros(F, sz)
    σm = zeros(F, sz)
    d = length(sz)

    ϵ1 = d == 2 ? ϵcore_ : ϵcore
    ϵ = device * (ϵ1) + (1 .- device) * ϵclad
    b = Zygote.Buffer(ϵ)
    copyto!(b, ϵ)
    a = model()
    ϵm = a * (ϵ1) + (1 .- a) * ϵclad
    place!(b, ϵm, round.(Int, (designs[1].bbox[1] - origin) / dx) .+ 1)
    ϵ = copy(b)

    if d == 3
        ϵ = sandwich(ϵ, round.(Int, [hbase, hwg, hclad] / dx)..., ϵbase, ϵclad)
    end

    p = apply(geometry_padding; ϵ, μ, σ, σm)
    global q = p
    p = apply(geometry_staggering, p)
end
# heatmap(q.σ.σEx) |> display
#=
Optimal design will maximize powers into port 1 and out of port 2. Monitor normals were set so both are positive. `metrics` function compute these figures of merit (FOM) quantities by a differentiable FDTD simulation . `loss` is then defined accordingly 
=#

function metrics(model, configs, dϵ=0; autodiff=true, history=nothing)
    p = make_geometry(model, configs, dϵ;)
    if !isnothing(history)
        ignore_derivatives() do
            push!(history, p[:ϵ])
        end
    end
    @unpack dx, dt, u0, field_padding, source_instances, monitor_instances = configs
    Enames = keys(u0.E)
    Hnames = keys(u0.H)

    # run simulation
    u = reduce((u, t) -> maxwell_update(u, p, t, dx, dt, field_padding, source_instances;), 0:dt:T[1], init=deepcopy(u0))

    u, fp, = reduce(T[1]+dt:dt:T[2], init=(u, 0, 0)) do (u, fp,), t
        fp_ = [
            begin
                E = [u[k, m] for k = Enames]
                H = [u[k, m] for k = Hnames]
                # fp = dict([ :E=>E, :H=>H])
                [E, H]

                # E' = E - (E ⋅ normal(m)) * E
                # H' = H - (H ⋅ normal(m)) * H
                # mp = (E' × H') ⋅ normal(m)
            end
            for m = monitor_instances
        ]
        # tp_ = map(fp_, monitor_instances) do (E, H), m
        #     mean(sum((E × H) .* normal(m))) * area(m)
        # end

        (
            maxwell_update(u, p, t, dx, dt, field_padding, source_instances),
            fp + dt / Δ[2] * fp_ * cispi(-2t),
            # tp + dt / Δ[2] * tp_,
        )
    end

    mp = map(fp, monitor_instances) do (E, H), m
        if polarization == :TE
            # global Ex, Hy
            Ex, Hy, Ez = invreframe(frame(m), vcat(E, H))
        else
        end
        ap, am = mode_decomp(mode, (; Ex, Hy))
        @show abs(ap)^2, abs(am)^2
        abs(ap)^2 + 0real(Ez[1])
        # real((Ex ⋅ nfp.Ex) ⋅ (Hy ⋅ nfp.Hy)) * nmp + 0real(Ez[1])
        # abs(Ex ⋅ nfp.Ex) ⋅ abs(Hy ⋅ nfp.Hy) * nmp + 0real(Ez[1])
    end
    # @show mp
    mp
end
# error()
y = getindex.(values(targets), :power)
function score(v,)
    # mae(F(1.2) * y, v / mp0)
    # mean(y - v / mp0)
    y[2] - v[2] / mp0
end

history = []
loss = model -> score(metrics(model, configs))
@show loss(model)
# iterations2d == 0 && error()
#=
We now do adjoint optimization. The first few iterations may show very little change but will pick up momentum
=#

opt = Adam(1)
opt_state = Flux.setup(opt, model)
# iterations2d = 66
# iterations2d = 400
# for i = 1:iterations2d
for i = 1:40
    println("$i")
    @time l, (dldm,) = withgradient(loss, model)
    l < 0 && break
    Flux.update!(opt_state, model, dldm)
    println(" $l\n")
end
# @save "$(@__DIR__)/2d_model_$(time()).bson" model

if dosave
    @save "$(@__DIR__)/$(name)_model.bson" model
    mask = model() .< 0.5
    @save "$(@__DIR__)/$(name)_mask.bson" mask
    Images.save("$(@__DIR__)/$(name).png", mask)
end
# error()

function runsave(model, configs; kw...)
    p = make_geometry(model, configs)
    @unpack u0, dx, dt, field_padding, source_instances, monitor_instances = configs
    @showtime global u = accumulate((u, t) ->
            maxwell_update!(deepcopy(u), p, t, dx, dt, field_padding, source_instances),
        0:dt:T[2], init=u0)

    # move to cpu for plotting
    if ongpu
        u, p, source_instances = cpu.((u, p, source_instances))
    end
    Hz = field.(u, :Hz)
    ϵEy = field(p, :ϵEy)
    dir = @__DIR__
    d = ndims(Hz[1])
    _name = "$(d)d_$name"
    # error()
    recordsim("$dir/$(_name).mp4", Hz, ;
        dt,
        field=:Hz,
        monitor_instances,
        source_instances,
        geometry=ϵEy,
        rel_lims=0.8,
        playback=1,
        axis1=(; title="$(replace( _name,"_"=>" ")|>titlecase)"),
        axis2=(; title="monitor powers"),
        kw...
    )

end
heatmap(model())

record = model -> runsave(model, configs)
record2d && record(model)
