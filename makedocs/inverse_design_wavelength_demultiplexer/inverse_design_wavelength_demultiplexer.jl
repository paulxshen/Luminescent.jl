using UnPack, LinearAlgebra, Random, StatsBase, Interpolations, Zygote, Jello, Flux
using Flux: mae, Adam
using Zygote: withgradient, Buffer
using BSON: @save, @load
using Optim: Options, minimizer
using Optim, CUDA
using GLMakie
# using Luminesce, LuminesceVisualization

dir = pwd()
include("$dir/src/main.jl")
include("$dir/../LuminesceVisualization.jl/src/main.jl")
include("$dir/scripts/startup.jl")

F = Float32
Random.seed!(1)

"training params"
dogpu = false
name = "inverse_design_wavelength_demultiplexer"
nbasis = 6 # complexity of design region

# loads design layout
@load "$(@__DIR__)/layout.bson" base signals ports designs
@load "$(@__DIR__)/modes.bson" λ1 λ2 modes1 modes2 lb ub dx hsub wwg hwg hclad ϵsub ϵclad ϵwg
λ = λ1 # wavelength [um]
f2 = λ1 / λ2
T1 = 1 + norm(ports[1].c - ports[2].c) / λ * sqrt(ϵwg) # simulation duration in [periods]
T2 = 1 / (f2 - 1)
ϵmin = ϵclad
hsub, wwg, hwg, hclad, dx, ub, lb = [hsub, wwg, hwg, hclad, dx, ub, lb] / λ

base = F.(base)
ϵsub, ϵwg, ϵclad = F.((ϵsub, ϵwg, ϵclad))
ϵdummy = sandwich(base, round.(Int, [hsub, hwg, hclad] / dx), [ϵsub, ϵwg, ϵclad])
sz0 = size(ϵdummy)

# "geometry generator model
# @load "model.bson" model
contrast = 10.0
model = Mask(round.(Int, designs[1].L / λ / dx) .+ 1, nbasis, contrast;)
model0 = deepcopy(model)


"monitors"
normal = [1, 0, 0] # normal 
δ = 0.1 / λ # margin
# (center, lower bound, upper bound; normal)
monitors = [Monitor([p.c / λ..., hsub], [0, -wwg / 2 - δ, -δ], [0, wwg / 2 + δ, hwg + δ]; normal=[p.n..., 0]) for p = ports]

# modal source
sources = []
c = signals[1].c / λ
c = [c..., hsub]
lb = [0, lb...]
ub = [0, ub...]
for (modes, f) = zip((modes1, modes2), (1, f2))
    @unpack Ex, Ey, Ez, = modes[1]
    Jy, Jz, Jx = map([Ex, Ey, Ez] / maximum(maximum.(abs, [Ex, Ey, Ez]))) do a
        reshape(a, 1, size(a)...)
    end
    push!(sources, Source(t -> cispi(2f * t), c, lb, ub; Jx, Jy, Jz))
end

boundaries = [] # unspecified boundaries default to PML
configs = setup(boundaries, sources, monitors, dx, sz0; ϵmin)
@unpack μ, σ, σm, dt, geometry_padding, geometry_splits, field_padding, source_instances, monitor_instances, u0, = configs
nt = round(Int, 1 / dt)

T = T1 + T2
t = 0:dt:T
v0 = zeros(F, length(monitor_instances))

t = F.(t)
dt, f2, T1, T2, T = F.((dt, f2, T1, T2, T))
if dogpu
    using CUDA, Flux
    @assert CUDA.functional()
    u0, model, base, μ, σ, σm, field_padding, source_instances =#, monitor_instances =
        gpu.((u0, model, base, μ, σ, σm, field_padding, source_instances))#, monitor_instances))
end

function make_geometry(model, base, μ, σ, σm)
    base_ = Buffer(base)
    base_[:, :] = base
    place!(base_, model(), round.(Int, designs[1].o / λ / dx) .+ 1)

    ϵ = sandwich(copy(base_), round.(Int, [hsub, hwg, hclad] / dx), [ϵsub, ϵwg, ϵclad])
    ϵ, μ, σ, σm = apply(geometry_padding; ϵ, μ, σ, σm)
    p = apply(geometry_splits; ϵ, μ, σ, σm)
end

function metrics(model; T1=T1, T2=T2, autodiff=true)
    p = make_geometry(model, base, μ, σ, σm)
    # run simulation
    _step = if autodiff
        step3
    else
        step3!
    end
    u = reduce((u, t) -> _step(u, p, t, dx, dt, field_padding, source_instances;), 0:dt:T1, init=deepcopy(u0))
    v = reduce(T1+dt:dt:T, init=(u, v0)) do (u, v), t
        _step(u, p, t, dx, dt, field_padding, source_instances),
        v + [1, cispi.(2 * F.([1, f2]))...] .* power.(monitor_instances, (u,))
    end[2]
    dt / F(T2) * abs.(v)

end
@show const tp = metrics(model, T1=1, T2=1, autodiff=false)[1] # total power
# error()

function loss(y)
    sum(abs, -y / tp)
end

p0 = make_geometry(model0, base, μ, σ, σm)
volume(cpu(p0[1][2]))
# error()

# gradient free optimization "Optim functions
x0, re = destructure(model)
f_ = m -> loss(metrics(m; autodiff=false))
# f_ = m -> loss(metrics(m;))
f = f_ ∘ re
x = deepcopy(x0)

CUDA.@allowscalar res = optimize(f, x,
    ParticleSwarm(; n_particles=10),
    Optim.Options(f_tol=0, iterations=0, show_every=1, show_trace=true))
xgf = minimizer(res)
x = deepcopy(xgf)
heatmap(cpu(re(x)()))
@show metrics(re(x))
@save "$(@__DIR__)/model.bson" model = re(x)
model = re(x)
# error()

# adjoint optimization
opt = Adam(0.2)
opt_state = Flux.setup(opt, model)
n = 1
for i = 1:n
    @time l, (dldm,) = withgradient(m -> loss(metrics(m)), model)
    Flux.update!(opt_state, model, dldm)
    println("$i $l")
end
heatmap(cpu(model()))
@save "$(@__DIR__)/model.bson" model
# error()

@show metrics(model)
function runsave(x)
    model = re(x)
    p = make_geometry(model, base, μ, σ, σm)
    @showtime u = accumulate((u, t) ->
            step3!(deepcopy(u), p, t, dx, dt, field_padding, source_instances),
        t, init=u0)

    # move to cpu for plotting
    global source_instances
    if dogpu
        u, p, source_instances = cpu.((u, p, source_instances))
    end

    f = :Ey
    Ey = field.(u, f)
    ϵEy = p[1][2]
    # error()
    recordsim("$(@__DIR__)/$(name).mp4", Ey, ;
        dt,
        field=f,
        monitor_instances,
        source_instances,
        geometry=ϵEy,
        elevation=60°,
        playback=1,
        axis1=(; title="$name\n$f"),
        axis2=(; title="monitor powers"),
    )

end

# runsave(x)