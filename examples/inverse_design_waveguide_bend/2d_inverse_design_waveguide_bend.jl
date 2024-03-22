using UnPack, LinearAlgebra, Random, StatsBase
using Zygote, Flux, Optim, CUDA, GLMakie, Jello
using Flux: mae, Adam
using Zygote: withgradient, Buffer
using BSON: @save, @load
using Optim: Options, minimizer
using AbbreviatedStackTraces
# using Jello, Luminescent, LuminescentVisualization
Random.seed!(1)

dir = pwd()
include("$dir/src/main.jl")
include("$dir/../LuminescentVisualization.jl/src/main.jl")
include("$dir/scripts/startup.jl")

# loads design layout
@load "$(@__DIR__)/layout.bson" base signals ports designs λ dx ϵsub ϵclad ϵcore

# training params"
iterations = 10
F = Float32
dogpu = false
name = "2d_inverse_design_signal_splitter"
nbasis = 5 # complexity of design region
contrast = 10
# rmin = round(Int, 0.1 / dx)
# rmin = :auto
rmin = nothing
init = nothing

T1 = Δ1 = 2 + 1.5norm(signals[1].c - ports[2].c) / λ * sqrt(ϵcore) # simulation duration in [periods] for signal to reach output ports
Δ2 = 4 # duration to record power at output ports
T2 = T = T1 + Δ2
ϵmin = ϵclad
dx, = [dx,] / λ

base = F.(base)
ϵsub, ϵcore, ϵclad = F.((ϵsub, ϵcore, ϵclad))
sz = size(base)

# "geometry generator model
# @load "$(@__DIR__)/model.bson" model
szd = Tuple(round.(Int, designs[1].L / λ / dx) .+ 1)
model = Blob(szd...; init, nbasis, contrast, rmin,)
# model.a .= [-1 -1 1 -1 -1; -1 1 1 -1 -1; 1 1 -1 -1 -1; -1 -1 -1 -1 -1; -1 -1 -1 -1 -1]
# model = randn(F, szd)
model0 = deepcopy(model)
heatmap(model())
# error()

# "boundaries"
boundaries = [] # unspecified boundaries default to PML

# "monitors"
normal = [1, 0, 0] # normal 
δ = 0.1 / λ # margin
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

μ = 1
σ = zeros(F, sz)
σm = zeros(F, sz)
configs = maxwell_setup(boundaries, sources, monitors, dx, sz; F, ϵmin)
@unpack dx, dt, geometry_padding, geometry_staggering, field_padding, source_instances, monitor_instances, u0, = configs

nt = round(Int, 1 / dt)
v0 = zeros(F, length(monitor_instances))

dx, dt, T1, T2, T = F.((dx, dt, T1, T2, T))

if dogpu
    using CUDA, Flux
    @assert CUDA.functional()
    u0, model, base, μ, σ, σm, field_padding, source_instances =
        gpu.((u0, model, base, μ, σ, σm, field_padding, source_instances))
end

function make_geometry(model, base, μ, σ, σm)
    base_ = Buffer(base)
    base_[:, :] = base
    # place!(base_, σ.(model), round.(Int, designs[1].o / λ / dx) .+ 1)
    place!(base_, model(), round.(Int, designs[1].o / λ / dx) .+ 1)
    mask = copy(base_)

    ϵ = mask * ϵcore + (1 .- mask) * ϵclad
    p = apply(geometry_padding; ϵ, μ, σ, σm)
    p = apply(geometry_staggering, p)
end

function metrics(model, ; T1=T1, T2=T2, autodiff=true)
    p = make_geometry(model, base, μ, σ, σm)
    # run simulation
    _step = if autodiff
        maxwell_update
    else
        maxwell_update!
    end
    u = reduce((u, t) -> _step(u, p, t, dx, dt, field_padding, source_instances;), 0:dt:T1, init=deepcopy(u0))
    v = reduce(T1+dt:dt:T2, init=(u, v0)) do (u, v), t
        _step(u, p, t, dx, dt, field_padding, source_instances),
        v + dt * power_flux.((u,), monitor_instances,)
    end[2] / Δ2
    abs.(v)
end
# @show const tp = metrics(model, T1=1, T2=2, autodiff=false)[1] # total power_flux
# error()

function loss(v)
    # sum(-v / tp)
    sum(-v)
end

p0 = make_geometry(model0, base, μ, σ, σm)

# adjoint optimization
opt = Adam(0.2)
opt_state = Flux.setup(opt, model)
# iterations = 66
# iterations = 400
for i = 1:iterations
    @time l, (dldm,) = withgradient(m -> loss(metrics(m)), model)
    Flux.update!(opt_state, model, dldm)
    println("$i $l")
end
@save "$(@__DIR__)/2d_model.bson" model
# error()

@show metrics(model)
function runsave(model)
    # model = re(x)
    p = make_geometry(model, base, μ, σ, σm)
    @showtime u = accumulate((u, t) ->
            maxwell_update!(deepcopy(u), p, t, dx, dt, field_padding, source_instances),
        0:dt:T2, init=u0)

    # move to cpu for plotting
    global source_instances
    if dogpu
        u, p, source_instances = cpu.((u, p, source_instances))
    end
    Ey = get.(u, :Ey)
    ϵEy = get(p, :ϵEy)
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