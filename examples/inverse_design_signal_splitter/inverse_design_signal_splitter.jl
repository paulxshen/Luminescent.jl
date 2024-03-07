using UnPack, LinearAlgebra, Random, StatsBase
using Zygote, Flux, Optim, CUDA, GLMakie, Jello
using Flux: mae, Adam
using Zygote: withgradient, Buffer
using BSON: @save, @load
using Optim: Options, minimizer
using AbbreviatedStackTraces
using Jello, Luminescent, LuminescentVisualization

# dir = pwd()
# include("$dir/src/main.jl")
# include("$dir/../LuminescentVisualization.jl/src/main.jl")
# include("$dir/scripts/startup.jl")

Random.seed!(1)

# training params"
F = Float32
dogpu = false
name = "inverse_design_signal_splitter"
nbasis = 6 # complexity of design region

# loads design layout
@load "$(@__DIR__)/layout.bson" base signals ports designs
@load "$(@__DIR__)/modes.bson" modes lb ub λ dx hsub wwg hwg hclad ϵsub ϵclad ϵwg
T1 = Δ1 = 1 + norm(ports[1].c - ports[2].c) / λ * sqrt(ϵwg) # simulation duration in [periods] for signal to reach output ports
Δ2 = 1 # duration to record power at output ports
T2 = T = T1 + Δ2
ϵmin = ϵclad
hsub, wwg, hwg, hclad, dx, ub, lb = [hsub, wwg, hwg, hclad, dx, ub, lb] / λ

base = F.(base)
ϵsub, ϵwg, ϵclad = F.((ϵsub, ϵwg, ϵclad))
ϵdummy = sandwich(base, round.(Int, [hsub, hwg, hclad] / dx), [ϵsub, ϵwg, ϵclad])
sz = size(ϵdummy)

# "geometry generator model
# @load "$(@__DIR__)/model.bson" model
contrast = 10.0
model = Mask(round.(Int, designs[1].L / λ / dx) .+ 1, nbasis, contrast; symmetries=2)
model0 = deepcopy(model)

# "boundaries"
boundaries = [] # unspecified boundaries default to PML

# "monitors"
normal = [1, 0, 0] # normal 
δ = 0.1 / λ # margin
monitors = [
    # (center, lower bound, upper bound; normal)
    Monitor([p.c / λ..., hsub], [0, -wwg / 2 - δ, -δ], [0, wwg / 2 + δ, hwg + δ]; normal=[p.n..., 0])
    for p = ports]

# modal source
@unpack Ex, Ey, Ez, = modes[1]
Jy, Jz, Jx = map([Ex, Ey, Ez] / maximum(maximum.(abs, [Ex, Ey, Ez]))) do a
    reshape(a, 1, size(a)...)
end
# GLMakie.volume(real(Jy))
c = [signals[1].c / λ..., hsub]
lb = [0, lb...]
ub = [0, ub...]
sources = [Source(t -> cispi(2t), c, lb, ub; Jx, Jy, Jz)]

μ = 1
σ = zeros(F, sz)
σm = zeros(F, sz)
configs = maxwell_setup(boundaries, sources, monitors, dx, sz; F, ϵmin)
@unpack dt, geometry_padding, geometry_staggering, field_padding, source_instances, monitor_instances, u0, = configs

nt = round(Int, 1 / dt)
v0 = zeros(F, length(monitor_instances))

dx, dt, T1, T2, T = F.((dx, dt, T1, T2, T))

if dogpu
    using CUDA, Flux
    @assert CUDA.functional()
    u0, model, base, μ, σ, σm, t, field_padding, source_instances =
        gpu.((u0, model, base, μ, σ, σm, t, field_padding, source_instances))
end

function make_geometry(model, base, μ, σ, σm)
    base_ = Buffer(base)
    base_[:, :] = base
    place!(base_, model(), round.(Int, designs[1].o / λ / dx) .+ 1)

    ϵ = sandwich(copy(base_), round.(Int, [hsub, hwg, hclad] / dx), [ϵsub, ϵwg, ϵclad])
    ϵ, μ, σ, σm = apply(geometry_padding; ϵ, μ, σ, σm)
    p = apply(geometry_staggering; ϵ, μ, σ, σm)
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
        v + dt * power_flux.(monitor_instances, (u,))
    end[2]
    abs.(v)
end
@show const tp = metrics(model, T1=1, T2=2, autodiff=false)[1] # total power_flux
# error()

function loss(v)
    sum(-v / tp)
end

p0 = make_geometry(model0, base, μ, σ, σm)
volume(cpu(p0[1][2]))
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
n = 60
for i = 1:n
    @time l, (dldm,) = withgradient(m -> loss(metrics(m)), model)
    Flux.update!(opt_state, model, dldm)
    println("$i $l")
end
heatmap(cpu(model()))
@save "$(@__DIR__)/model.bson" model
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
    Ey = map(u) do u
        u[1][2]
    end
    ϵEy = p[1][2]
    dir = @__DIR__
    # error()
    recordsim("$dir/$(name).mp4", Ey, ;
        dt,
        field=:Ey,
        monitor_instances,
        source_instances,
        geometry=ϵEy,
        elevation=60°,
        playback=1,
        axis1=(; title="$(replace( name,"_"=>" ")|>titlecase)"),
        axis2=(; title="monitor powers"),
    )

end

runsave(model)