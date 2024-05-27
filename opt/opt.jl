#=
We do inverse design of a compact photonic waveguide bend to demonstrate workflow of FDTD adjoint optimization. First, we seed the design using 2d TE adjoint simulations which serve as fast approximations. Optionlly, we finetune the resulting design in full blown 3d adjoint simulations.

=#

using UnPack, LinearAlgebra, Random, StatsBase, Dates, Colors, DataStructures
using Zygote, Flux, CUDA, GLMakie, Jello, ImageTransformations
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
# dir = @__DIR__
dir = pwd()
@load "$(dir)/_prob.bson" name signals ports designs λc dx masks path_length targets margins lmin wavelengths modes
#  ϵbase ϵclad ϵcore hbase hwg hclad
signals, ports, targets = (signals, ports, targets) .|> SortedDict
polarization = :TE

# error()
for s = [:device, :guess]
    # for name = ["device", "guess"]
    r = Rsvg.handle_new_from_data(masks[s][:svg])
    # r = Rsvg.handle_new_from_file("$(@__DIR__)/$name.svg")
    d = Rsvg.handle_get_dimensions(r)
    cs = Cairo.CairoImageSurface(1 * (d.width), 1 * (d.height), Cairo.FORMAT_ARGB32)
    c = Cairo.CairoContext(cs)
    Rsvg.handle_render_cairo(c, r)
    Cairo.write_to_png(cs, "$(@__DIR__)/$s.png")
end

masks = Dict([k => imresize(F.(convert.(Gray, FileIO.load("$(@__DIR__)/$k.png") |> transpose),), round.(Int, (masks[k].bbox[2] - masks[k].bbox[1]) / dx) |> Tuple) .> tol for k in [:device, :guess]])
@unpack device, guess = masks
# heatmap(masks.device) |> display
lm, rm = round.(margins / dx)
device = pad(device, 0, lm, rm)
origin = masks.device.bbox[1] - lm * dx
#=
We initialize a Jello.jl Blob object which will generate geometry of design region. Its parameters will get optimized during adjoint optimization. We initialize it with a straight slab connecting input to output port.
=#

bbox = designs[1].bbox
L = bbox[2] - bbox[1]
szd = Tuple(round.(Int, L / dx)) # design region size
model = Blob(szd;
    init=guess,
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
Δ[1] = 2 + path_length / λc * sqrt(ϵcore) # simulation duration in [periods] for signal to reach output ports
Δ[2] = 2 # duration to record power at output ports
T = cumsum(Δ)

#=
We set boundary conditions, sources , and monitor. The modal source profile is obtained from external mode solver , in our case VectorModesolver.jl . Please refer to guide section of docs website for details . To get an approximate  line source for use in 2d from the cross section profile , we sum and collapse it along its height axis
=#

boundaries = [] # unspecified boundaries default to PML
sig = (signals |> values |> first)
wm = sig.size[1]
monitors = [
    # (center, lower bound, upper bound; normal)
    begin
        c = (p.center - origin) / λc
        n = p.normal
        if k in keys(signals)
            c -= margin_offset * n
            c -= source_monitor_offset * n
            n *= -1
        else
            c -= source_monitor_offset * n
        end
        L = (p.endpoints[2] - p.endpoints[1])
        L = L / norm(L) * wm / λc
        Monitor(c, L, n)
    end
    for (k, p) = ports |> pairs
]
# monitors = [monitors[1]]

device, guess, T, Δ, ϵbase, ϵcore, ϵclad, dx, λc =
    F.((device, guess, T, Δ, ϵbase, ϵcore, ϵclad, dx, λc))

# modal source
sources = []
_modes = []
for (port, sig) = signals |> pairs
    @unpack center, endpoints, wavelengths, mode_numbers = sig
    for (λ, mn) in zip(wavelengths, mode_numbers)
        v = (
            findfirst(mode_info) do v
                isapprox(λ, v.wavelength) && port ⊂ v.ports
            end
        )
        mode = v.mode[mn+1]
        sz = v.size
        eps = v.eps

        n = sig.normal
        c = (sig.center - origin) / λc + margin_offset * sig.normal
        L = (sig.endpoints[2] - sig.endpoints[1]) / λc
        mode = (; [k => transpose(complex.(stack.(v)...)) |> F for (k, v) in mode |> pairs]...)
        ϵ = eps |> stack |> transpose .|> F

        sz = round(sig.size / dx) |> Tuple
        mode = (; Pair.(keys(mode), [resize(v, size(v) - 1) for v in values(mode)])...)
        ϵ = resize(ϵ, sz)

        mode, ϵmode = collapse_mode(mode, ϵ)
        mode = normalize_mode(mode, dx / λc)
        i = round(Int, size(ϵmode, 1) / 2)
        global ϵcore_ = ϵmode[i]
        ϵslice = maximum(ϵ, dims=2) |> vec
        ϵslice = min.(ϵslice, ϵcore_)
        mode0 = deepcopy(mode)

        mode, a, b = calibrate_mode(mode, ϵslice, dx / λc; name="1")
        # error()
        mode = normalize_mode(mode, dx / λc)
        mode_, mp0, = calibrate_mode(mode, ϵslice, dx / λc; name="2")
        mode /= sqrt(mp0)
        @show mp0
        # @unpack Ex, Ez = mode

        # error()
        # sources = [ModalSource(t -> cispi(2t), c, L, n; Jx=mode.Ex, Jz=mode.Ez)]
        push!(sources, ModalSource(t -> cispi(2t * λc / λ), mode, c, L, n))
    end
end

ϵmin = ϵclad
sz = size(device)

d = 2
prob = maxwell_setup(boundaries, sources, monitors, dx / λc, sz; F, ϵmin, polarization, wavelengths)
@unpack dt, sz, geometry_padding, subpixel_averaging, field_padding, source_instances, monitor_instances, u0, polarization = prob
nt = round(Int, 1 / dt)
A = area.(monitor_instances)
# nmp = nfp.Ex ⋅ nfp.Hy * A[1] / length(Ex) / 2

if ongpu
    using Flux
    # using CUDA
    # @assert CUDA.functional()
    u0, model, device, μ, σ, σm, field_padding, source_instances =
        gpu.((u0, model, device, μ, σ, σm, field_padding, source_instances))
    merge!(prob, (; u0, field_padding, source_instances))
end

#=
We define a geometry update function that'll be called each adjoint iteration. It calls geometry generator model to generate design region which gets placed onto mask of device features.
    =#

ϵ2 = d == 2 ? ϵcore_ : ϵcore
ϵ1 = ϵclad
model_origin = round.(Int, (designs[1].bbox[1] - origin) / dx) .+ 1
function make_geometry(mask, ϵ1, ϵ2, model, origin)#; make3d=false)
    sz = size(mask)
    μ = ones(F, sz)
    σ = zeros(F, sz)
    σm = zeros(F, sz)
    d = length(sz)

    ϵ = mask * (ϵ2) + (1 .- mask) * ϵ1
    b = Zygote.Buffer(ϵ)
    copyto!(b, ϵ)
    a = model()
    ϵm = a * (ϵ2) + (1 .- a) * ϵ1
    place!(b, ϵm, origin)
    ϵ = copy(b)

    if d == 3
        ϵ = sandwich(ϵ, round.(Int, [hbase, hwg, hclad] / dx)..., ϵbase, ϵclad)
    end

    (; ϵ, μ, σ, σm)
end
# heatmap(q.σ.σEx) |> display
#=
Optimal design will maximize powers into port 1 and out of port 2. Monitor normals were set so both are positive. `metrics` function compute these figures of merit (FOM) quantities by a differentiable FDTD simulation . `loss` is then defined accordingly 
=#

y = getindex.(values(targets), :power)
ij = [(t.port[2] |> Int, findfirst(isapprox(t.wavelength), wavelengths)) for t in values(targets)]
function score(a)
    mae([abs(a[1, i, j])^2 for (i, j,) = ij], y)
end

function loss(model)
    prob[:geometry] = make_geometry(device, ϵ1, ϵ2, model, model_origin)
    score(solve(prob))
end
@show loss(model)
# iterations2d == 0 && error()
#=
We now do adjoint optimization. The first few iterations may show very little change but will pick up momentum
=#

opt = Adam(1)
opt_state = Flux.setup(opt, model)
# iterations2d = 66
# iterations2d = 400
for i = 1:iterations2d
    # for i = 1:40
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

function runsave(model, prob; kw...)
    p = make_geometry(model, prob)
    @unpack u0, dx, dt, field_padding, source_instances, monitor_instances = prob
    @showtime global u = accumulate((u, t) ->
            maxwell_update!(deepcopy(u), p, t, dx, dt, field_padding, source_instances),
        0:dt:T[2], init=u0)

    # move to cpu for plotting
    if ongpu
        u, p, source_instances = cpu.((u, p, source_instances))
    end
    Hz = field.(u, :Hz)
    ϵEy = field(p, :ϵyy)
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

record = model -> runsave(model, prob)
record2d && record(model)
