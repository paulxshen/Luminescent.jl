"""
We do inverse design of a compact photonic waveguide bend to demonstrate workflow of FDTD adjoint optimization. First, we seed the design in 2d TE adjoint simulations which fast approximate 3d propagation. Finally, we finetune the design in full blown 3d adjoint simulations.
"""

using UnPack, LinearAlgebra, Random, StatsBase
using Zygote, Flux, Optim, CUDA, GLMakie, Jello
using Flux: mae, Adam
using Zygote: withgradient, Buffer
using BSON: @save, @load
using Optim: Options, minimizer
using AbbreviatedStackTraces
# using Jello, Luminescent, LuminescentVisualization
Random.seed!(1)

# if running directly without module
dir = pwd()
include("$dir/src/main.jl")
include("$dir/../LuminescentVisualization.jl/src/main.jl")
include("$dir/scripts/startup.jl")

"""
We load design layout which includes a 2d static_mask of static waveguide geometry as well as variables with locations of ports, signals, design regions and material properties.
"""

@load "$(@__DIR__)/layout.bson" static_mask signals ports designs λ dx ϵsub ϵclad ϵcore hsub hwg hclad

# training params"
iterations = 10
F = Float32
dogpu = false
name = "2d_inverse_design_signal_splitter"

"""
We initialize a Jello.jl Blob object which will generate geometry of design region. Its parameters will get optimized during adjoint optimization. 
"""

nbasis = 5 # complexity of design region
contrast = 10 # edge sharpness 
rmin = nothing
init = 1 # uniform slab
szd = Tuple(round.(Int, designs[1].L / λ / dx) .+ 1)
model = Blob(szd...; init, nbasis, contrast, rmin,)
model0 = deepcopy(model)
heatmap(model())

"""
We set key time intervals. The signal must first propagate to port 2 after which all port power fluxes will get monitored
"""
Δ = zeros(2)
Δ[1] = 2 + 1.6norm(signals[1].c - ports[2].c) / λ * sqrt(ϵcore) # simulation duration in [periods] for signal to reach output ports
Δ[2] = 2 # duration to record power at output ports
T = cumsum(Δ)


# "geometry generator model
# @load "$(@__DIR__)/model.bson" model
# model.a .= [-1 -1 1 -1 -1; -1 1 1 -1 -1; 1 1 -1 -1 -1; -1 -1 -1 -1 -1; -1 -1 -1 -1 -1]
# model = randn(F, szd)
# error()

# "boundaries"
boundaries = [] # unspecified boundaries default to PML
monitors = [
    # (center, lower bound, upper bound; normal)
    Monitor(p.c / λ, p.lb / λ, p.ub / λ; normal=p.n)
    for p = ports
]

# modal source
@unpack Ex, Ey, Ez, = signals[1].modes[1]
Jy, Jx = map([Ex, Ez]) do a
    transpose(sum(a, dims=2))
end
Jy, Jx = [Jy, Jx] / maximum(maximum.(abs, [Jy, Jx]))
c = signals[1].c / λ
lb_ = [0, signals[1].lb[1]]
ub_ = [0, signals[1].ub[1]]
sources = [Source(t -> cispi(2t), c, lb_, ub_; Jx, Jy,)]

ϵmin = ϵclad
dx, = [dx,] / λ

static_mask = F.(static_mask)
ϵsub, ϵcore, ϵclad = F.((ϵsub, ϵcore, ϵclad))
sz = size(static_mask)

configs = maxwell_setup(boundaries, sources, monitors, dx, sz; F, ϵmin)
@unpack dx, dt, sz, geometry_padding, geometry_staggering, field_padding, source_instances, monitor_instances, u0, = configs

nt = round(Int, 1 / dt)
port_powers0 = zeros(F, length(monitor_instances))


if dogpu
    using CUDA, Flux
    @assert CUDA.functional()
    u0, model, static_mask, μ, σ, σm, field_padding, source_instances =
        gpu.((u0, model, static_mask, μ, σ, σm, field_padding, source_instances))
end

"""
We define a geometry update function that'll be called each adjoint iteration. It calls geometry generator model to generate design region which gets placed onto mask of static features.
    """
function make_geometry(model, static_mask, configs; make3d=false)
    @unpack sz, geometry_padding, geometry_staggering = configs
    μ = 1
    σ = zeros(F, sz)
    σm = zeros(F, sz)
    base_ = Buffer(static_mask)
    base_[:, :] = static_mask
    # place!(base_, σ.(model), round.(Int, designs[1].o / λ / dx) .+ 1)
    place!(base_, model(), round.(Int, designs[1].o / λ / dx) .+ 1)
    mask = copy(base_)
    ϵ = mask * ϵcore + (1 .- mask) * ϵclad

    if make3d
        ϵ = sandwich(ϵ, round.(Int, [hsub, hwg, hclad] / dx)..., ϵsub, ϵclad)
    end

    p = apply(geometry_padding; ϵ, μ, σ, σm)
    p = apply(geometry_staggering, p)
end

"""
Optimal design will maximize powers into port 1 and out of port 2. Monitor normals were set so both are positive. `metrics` function compute these figures of merit (FOM) quantities by a differentiable FDTD simulation . `loss` is then defined accordingly 
"""
function metrics(p, configs; autodiff=true)
    @unpack u0, field_padding, source_instances, monitor_instances = configs
    # run simulation
    _step = if autodiff
        maxwell_update
    else
        maxwell_update!
    end
    u = reduce((u, t) -> _step(u, p, t, dx, dt, field_padding, source_instances;), 0:dt:T[1], init=deepcopy(u0))
    port_powers = reduce(T[1]+dt:dt:T[2], init=(u, port_powers0)) do (u, port_powers), t
        _step(u, p, t, dx, dt, field_padding, source_instances),
        port_powers + dt * power_flux.((u,), monitor_instances,)
    end[2] / Δ[2]
    abs.(port_powers)
end
# @show const tp = metrics(model, T[1]=1, T[2]=2, autodiff=false)[1] # total power_flux
# error()

function score(port_powers)
    sum(-port_powers)
end

# p0 = make_geometry(model0, static_mask, μ, σ, σm)
loss = model -> score(metrics(make_geometry(model, static_mask, configs), configs))

"""
We now do adjoint optimization. The first few iterations may show very little change but will pick up momentum
"""
opt = Adam(0.5)
opt_state = Flux.setup(opt, model)
# iterations = 66
# iterations = 400
for i = 1:iterations
    @time l, (dldm,) = withgradient(loss, model)
    Flux.update!(opt_state, model, dldm)
    println("$i $l")
end
@save "$(@__DIR__)/2d_model_$(now()).bson" model
# error()

"""
We do a simulation movie using optimized geometry
"""
@show metrics(model)
function runsave(model)
    p = make_geometry(model, static_mask, μ, σ, σm)
    @showtime u = accumulate((u, t) ->
            maxwell_update!(deepcopy(u), p, t, dx, dt, field_padding, source_instances),
        0:dt:T[2], init=u0)

    # move to cpu for plotting
    global source_instances
    if dogpu
        u, p, source_instances = cpu.((u, p, source_instances))
    end
    Ey = field.(u, :Ey)
    ϵEy = field(p, :ϵEy)
    dir = @__DIR__
    # error()
    recordsim("$dir/$(name).mp4", Ey, ;
        dt,
        field=:Ey,
        monitor_instances,
        source_instances,
        geometry=ϵEy,
        rel_lims=0.5,
        playback=1,
        axis1=(; title="$(replace( name,"_"=>" ")|>titlecase)"),
        axis2=(; title="monitor powers"),
    )

end

# heatmap(cpu(model()))
runsave(model)

"""
We now finetune our design in 3d by starting off with optimized model from 2d. We make 3d geometry simply by sandwiching thickened 2d mask between lower substrate and upper clad layers. 
"""
ϵdummy = sandwich(static_mask, round.(Int, [hsub, hwg, hclad] / dx)..., ϵsub, ϵclad)
sz = size(ϵdummy)
model2d = deepcopy(model)


# "monitors"
δ = 0.1 # margin
monitors = [Monitor([p.c / λ..., hsub / λ], [p.lb / λ..., -δ / λ], [p.ub / λ..., hwg / λ + δ / λ]; normal=[p.n..., 0]) for p = ports]

# modal source
@unpack Ex, Ey, Ez, = signals[1].modes[1]
Jy, Jz, Jx = map([Ex, Ey, Ez] / maximum(maximum.(abs, [Ex, Ey, Ez]))) do a
    reshape(a, 1, size(a)...)
end
c = [signals[1].c / λ..., hsub / λ]
lb = [0, signals[1].lb...] / λ
ub = [0, signals[1].ub...] / λ
sources = [Source(t -> cispi(2t), c, lb, ub; Jx, Jy, Jz)]

configs = maxwell_setup(boundaries, sources, monitors, dx, sz; F, ϵmin)

loss = model -> score(metrics(make_geometry(model, static_mask, configs; make3d=true), configs))