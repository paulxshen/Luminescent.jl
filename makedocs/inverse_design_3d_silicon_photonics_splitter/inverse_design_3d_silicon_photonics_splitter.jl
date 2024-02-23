using UnPack, LinearAlgebra, Random, StatsBase, Interpolations, Zygote, Jello, Flux
using Flux: mae, Adam
using Zygote: withgradient
using BSON: @save, @load
using Optim: Options, minimizer
using Optim

using GLMakie
dir = pwd()
include("$dir/src/main.jl")
include("$dir/../FDTDToolkit.jl/src/main.jl")
include("$dir/scripts/startup.jl")

# dir = "FDTDEngine.jl"
# using FDTDEngine,FDTDToolkit

F = Float32
dogpu = true
# dogpu = false
if dogpu# &&
    using CUDA, Flux
    @assert CUDA.functional()
end
Random.seed!(1)

"training params"
name = "silicon_photonics_splitter"
nx = 24
T = 18.0f0 # simulation duration in [periods]
nbasis = 4 # complexity of design region
ϵmin = 2.26
λ = 1.55f0 # wavelength [um]

# loads design layout
@load "$(@__DIR__)/layout.bson" layout
@unpack wwg, base, sources, ports, designs, dx, = layout
base = F.(base)
hcore = 0.22f0 / λ
hclad = hsub = 0.2f0 / λ
ϵsub = ϵclad = 2.25f0
ϵcore = 12.25f0
ϵdummy = sandwich(base, round.(Int, [hsub, hcore, hclad] / dx), [ϵsub, ϵcore, ϵclad])
sz0 = size(ϵdummy)

"geometry generator model"
contrast = 10.0f0
model = Mask(round.(Int, designs[1].L / dx) .+ 1, nbasis, contrast)
model0 = deepcopy(model)

"boundaries"
boundaries = [] # unspecified boundaries default to PML

"monitors"
z = hsub + hcore / 2
monitors = vcat(
    [Monitor([p.c..., z], [0, wwg, hcore], [p.n..., 0]) for p = ports],
    monitors_on_box([designs[1].c..., z], [designs[1].L..., hcore])
)


# modal source
@load "$(@__DIR__)/mode.bson" mode
@unpack Ex, Ey, Ez, bounds = mode
bounds = f32(bounds / λ)
lb = getindex.(bounds, 1)
ub = getindex.(bounds, 2)
_dx = mode.dx / λ
Ex, Ey, Ez = [Ex, Ey, Ez] / maximum(maximum.(abs, [Ex, Ey, Ez]))
sy, sz = size.((Ex,), (1, 2))
y = range(F.(bounds[1])..., sy)
z = range(F.(bounds[2])..., sz)
gEx = LinearInterpolation((y, z), F.(Ez))
gEy = LinearInterpolation((y, z), F.(Ex))
gEz = LinearInterpolation((y, z), F.(Ey))

c = [sources[1].c..., hsub]
lb = [0, lb...]
ub = [0, ub...]
sources = [
    Source(t -> cos(F(2π) * t), c, lb, ub;
        Jx=(x, y, z) -> gEx(y, z),
        Jy=(x, y, z) -> gEy(y, z),
        Jz=(x, y, z) -> gEz(y, z),
    ),
]

configs = setup(boundaries, sources, monitors, dx, sz0; ϵmin, T)
@unpack μ, σ, σm, dt, geometry_padding, geometry_splits, field_padding, source_instances, monitor_instances, u0, = configs
nt = round(Int, 1 / dt)

t = 0:dt:T
y0 = zeros(F, length(monitor_instances))
if dogpu# &&
    u0, model, base, μ, σ, σm, t, field_padding, source_instances =#, monitor_instances =
        gpu.((u0, model, base, μ, σ, σm, t, field_padding, source_instances))#, monitor_instances))
end

function make_geometry(model, base, μ, σ, σm)
    base_ = Buffer(base)
    base_[:, :] = base
    place!(base_, model(), round.(Int, designs[1].o / dx) .+ 1)

    ϵ = sandwich(copy(base_), round.(Int, [hsub, hcore, hclad] / dx), [ϵsub, ϵcore, ϵclad])
    ϵ, μ, σ, σm = apply(geometry_padding; ϵ, μ, σ, σm)
    p = apply(geometry_splits; ϵ, μ, σ, σm)
end

function metrics(model;)
    p = make_geometry(model, base, μ, σ, σm)
    # run simulation
    u = reduce((u, t) -> step(u, p, t, dx, dt, field_padding, source_instances;), 0:dt:T-1, init=u0)
    y = reduce(
        ((u, y), t) -> (
            step(u, p, t, dx, dt, field_padding, source_instances),
            y + dt * power.(monitor_instances, (u,))),
        T-1+dt:dt:T,
        init=(u, y0))[2]

end
function loss(y)
    dot([0, -1, -1, 1, -1, 1, 1, 1, 1], y) +
    abs(y[2] - y[3])
end

tp = abs(metrics(model)[1]) # total power
p0 = make_geometry(model0, base, μ, σ, σm)
# volume(p0[1][3])

"Optim functions"
x0, re = destructure(model)
f_ = m -> loss(metrics(m))
f = f_ ∘ re
x = deepcopy(x0)
@show f(x)


opt = Adam(0.2)
m = re(x)
opt_state = Flux.setup(opt, m)
n = 12
for i = 1:n
    @time l, (dldm,) = withgradient(m -> loss(metrics(m)), m)
    Flux.update!(opt_state, m, dldm)
    println("$i $l")
end
# x_adjoint = re(minimizer(res))

# move to cpu for plotting
if dogpu
    u0, model, base, μ, σ, σm, t, field_padding, source_instances =#, monitor_instances =
        cpu.((u0, model, base, μ, σ, σm, t, field_padding, source_instances))#, monitor_instances))
end

p = make_geometry(m, base, μ, σ, σm)
volume(cpu(p[1][3]))

function runsave(x)
    m = re(x)
    p = make_geometry(m, base, μ, σ, σm)
    @showtime u = accumulate((u, t) ->
            step!(deepcopy(u), p, t, dx, dt, field_padding, source_instances),
        t, init=u0)

    global source_instances
    if dogpu
        u, p, source_instances = cpu.((u, p, source_instances))
    end
    Ez = map(u) do u
        u[1][3]
    end
    ϵz = p[1][3]
    dir = @__DIR__
    # error()
    recordsim("$dir/$(name)_nres_$nx.mp4", collect.(Ez), ;
        dt,
        monitor_instances,
        source_instances,
        geometry=ϵz,
        elevation=30°,
        playback=1,
        axis1=(; title="$name\nEz"),
        axis2=(; title="monitor powers"),
    )

end

runsave(x)