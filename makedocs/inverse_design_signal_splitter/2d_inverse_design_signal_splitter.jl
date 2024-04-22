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
@load "$(@__DIR__)/layout.bson" mask signals ports designs
@load "$(@__DIR__)/modes.bson" modes lb ub λ dx hsub wwg hwg hclad ϵbase ϵclad ϵcore

# training params"
F = Float32
dogpu = false
name = "inverse_design_signal_splitter"
nbasis = 10 # complexity of design region
contrast = 10
# rmin = round(Int, 0.1 / dx)
# rmin = :auto
rmin = nothing
init = nothing

T[1] = Δ[1] = 1 + 1.3norm(ports[1].c - ports[2].c) / λ * sqrt(ϵcore) # simulation duration in [periods] for signal to reach output ports
Δ[2] = 1 # duration to record power at output ports
T[2] = T = T[1] + Δ[2]
ϵmin = ϵclad
hsub, wwg, hwg, hclad, dx, ub, lb = [hsub, wwg, hwg, hclad, dx, ub, lb] / λ

mask = F.(mask)
ϵbase, ϵcore, ϵclad = F.((ϵbase, ϵcore, ϵclad))
sz = size(mask)

# "geometry generator model
# @load "$(@__DIR__)/model.bson" model
model = RealBlob((round.(Int, designs[1].L / λ / dx) .+ 1)...;
    init, nbasis, contrast, rmin, symmetry_dims=2)
model0 = deepcopy(model)

# "boundaries"
boundaries = [] # unspecified boundaries default to PML

# "monitors"
normal = [1, 0, 0] # normal 
δ = 0.1 / λ # margin
monitors = [
    # (center, lower bound, upper bound; normal)
    Monitor(p.c / λ, [0, -wwg / 2 - δ], [0, wwg / 2 + δ]; normal=p.n)
    for p = ports]

# modal source
@unpack Ex, Ey, Ez, = modes[1]
Jy, Jx = map([Ex, Ez]) do a
    transpose(sum(a, dims=2))
end
Jy, Jx = [Jy, Jx] / maximum(maximum.(abs, [Jy, Jx]))
c = signals[1].c / λ
lb_ = [0, lb[1]]
ub_ = [0, ub[1]]
sources = [Source(t -> cispi(2t), c, lb_, ub_; Jx, Jy,)]

μ = 1
σ = zeros(F, sz)
σm = zeros(F, sz)
configs = maxwell_setup(boundaries, sources, monitors, dx, sz; F, ϵmin)
@unpack dx, dt, geometry_padding, geometry_staggering, field_padding, source_instances, monitor_instances, u0, = configs

nt = round(Int, 1 / dt)
port_fluxes0 = zeros(F, length(monitor_instances))

dx, dt, T[1], T[2], T = F.((dx, dt, T[1], T[2], T))

if dogpu
    using CUDA, Flux
    @assert CUDA.functional()
    u0, model, mask, μ, σ, σm, field_padding, source_instances =
        gpu.((u0, model, mask, μ, σ, σm, field_padding, source_instances))
end

function make_geometry(model, mask, μ, σ, σm)
    base_ = Buffer(mask)
    base_[:, :] = mask
    place!(base_, model(), round.(Int, designs[1].o / λ / dx) .+ 1)
    mask = copy(base_)

    ϵ = mask * ϵcore + (1 .- mask) * ϵclad
    p = apply(geometry_padding; ϵ, μ, σ, σm)
    p = apply(geometry_staggering, p)
end

function metrics(model, ; T[1]=T[1], T[2]=T[2], autodiff=true)
    p = make_geometry(model, mask, μ, σ, σm)
    # run simulation
    _step = if autodiff
        maxwell_update
    else
        maxwell_update!
    end
    u = reduce((u, t) -> _step(u, p, t, dx, dt, field_padding, source_instances;), 0:dt:T[1], init=deepcopy(u0))
    v = reduce(T[1]+dt:dt:T[2], init=(u, port_fluxes0)) do (u, v), t
        _step(u, p, t, dx, dt, field_padding, source_instances),
        v + dt * power.((u,), monitor_instances,)
    end[2]
    abs.(v)
end
# @show const tp = metrics(model, T[1]=1, T[2]=2, autodiff=false)[1] # total power
# error()

function loss(v)
    # sum(-v / tp)
    sum(-v)
end

p0 = make_geometry(model0, mask, μ, σ, σm)
# volume(cpu(p0[1][2]))
# error()

# gradient free optimization "Optim functions
# x0, re = destructure(model)
# f_ = m -> loss(metrics(m; autodiff=false))
# f = f_ ∘ re
# x = deepcopy(x0)

# CUDA.@allowscalar res = optimize(f, x,
#     ParticleSwarm(; n_particles=10),
#     Optim.Options(f_tol=0, iterations=10, show_every=1, show_trace=true))
# xgf = minimizer(res)
# x = deepcopy(xgf)
# heatmap(cpu(re(x)()))
# @show metrics(re(x))
# @save "$(@__DIR__)/model.bson" model = re(x)
# model = re(x)
# error()

# adjoint optimization
opt = Adam(0.2)
opt_state = Flux.setup(opt, model)
n = 66
n = 120
n = 400
for i = 1:n
    @time l, (dldm,) = withgradient(m -> loss(metrics(m)), model)
    Flux.update!(opt_state, model, dldm)
    println("$i $l")
end
@save "$(@__DIR__)/2d_model.bson" model
# error()

@show metrics(model)
function runsave(model)
    # model = re(x)
    p = make_geometry(model, mask, μ, σ, σm)
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
        playback=1,
        axis1=(; title="$(replace( name,"_"=>" ")|>titlecase)"),
        axis2=(; title="monitor powers"),
    )

end

# heatmap(cpu(model()))
runsave(model)