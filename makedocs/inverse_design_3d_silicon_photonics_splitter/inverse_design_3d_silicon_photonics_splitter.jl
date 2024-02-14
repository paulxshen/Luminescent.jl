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
Random.seed!(1)

"training params"
name = "silicon_photonics_splitter"
nres = 16
T = 18.0f0 # simulation duration in [periods]
nbasis = 4 # complexity of design region
ϵmin = 1.5
λ = 1.55f0 # wavelength [um]

# loads design layout
@load "$(@__DIR__)/layout.bson" layout
@unpack wwg, base, sources, ports, designs, dx = layout
hcore = 0.22f0 / λ
hclad = hsub = 0.2f0 / λ
ϵsub = ϵclad = 2.25f0
ϵcore = 12.25f0
ϵdummy = sandwich(base, round.(Int, [hsub, hcore, hclad] / dx), [ϵsub, ϵcore, ϵclad])
sz0 = size(ϵdummy)

"geometry generator model"
contrast = 10.0f0
nbasis = 4
model = Mask(round.(Int, designs[1].L / dx), nbasis, contrast)
model0 = deepcopy(model)

"boundaries"
boundaries = [] # unspecified boundaries default to PML

"monitors"
z = [hsub, hsub + hcore]
monitors = vcat(
    [Monitor([xy..., z], n) for (xy, n) = ports],
    monitors_on_box([designs[1].c..., hsub], designs[1].L)
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

function make_geometry(model, base, μ, σ, σm)
    base = place!(F.(base), model(), round.(Int, c / dx) .+ 1)
    ϵ = sandwich(base, round.(Int, [hsub, hcore, hclad] / dx), [ϵsub, ϵcore, ϵclad])
    ϵ, μ, σ, σm = apply(geometry_padding; ϵ, μ, σ, σm)
    p = apply(geometry_splits; ϵ, μ, σ, σm)
end
p0 = make_geometry(model0, base, μ, σ, σm)
# volume(p0[1][3])

function metrics(model, T=T;)
    p = make_geometry(model, base, μ, σ, σm)
    # run simulation
    u = reduce((u, t) -> step!(u, p, t, dx, dt, field_padding, source_instances;), 0:dt:T-1, init=deepcopy(u0))
    y = reduce(((u, y), t) -> (
            step!(u, p, t, dx, dt, field_padding, source_instances),
            y + dt * power.(monitor_instances, (u,))
        ),
        T-1+dt:dt:T,
        init=(u, zeros(F, length(monitor_instances))))[2]
end
tp = abs(metrics(model, 2)[1]) # total power

function loss(y)
    dot([0, -1, -1, 1, -1, 1, 1, 1, 1], y) +
    abs(y[2] - y[3])
end

"Optim functions"
x0, re = destructure(model)
f_ = loss ∘ metrics
f = f_ ∘ re
x = deepcopy(x0)
@show f(x)


opt = Adam(0.2)
m = re(x)
opt_state = Flux.setup(opt, m)
n = 2
for i = 1:n
    @time l, (dldm,) = withgradient(m -> loss(metrics(m;)), m)
    Flux.update!(opt_state, m, dldm)
    println("$i $l")
end
# x_adjoint = re(minimizer(res))

function runsave(x)
    model = re(x)
    p = make_geometry(model, base, μ, σ, σm)
    t = 0:dt:T
    sol = similar([u0], length(t))
    @showtime sol = accumulate!((u, t) -> step!(u, p, t, dx, dt, field_padding, source_instances), sol, t, init=u0)

    # make movie
    Ez = map(sol) do u
        u[1][3]
    end
    ϵz = p[1][3]
    dir = @__DIR__
    ° = π / 180
    recordsim(Ez, ϵz, configs, "$dir/$(name)_nres_$nres.mp4", title="$name"; elevation=60°, playback=1, bipolar=true)

end

runsave(xgf)