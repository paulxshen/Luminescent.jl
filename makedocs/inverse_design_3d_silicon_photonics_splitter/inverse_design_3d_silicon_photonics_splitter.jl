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
Courant = 1.5 * 0.7 / √3# Courant number
C = 1000 # scale loss
λ = 1.55f0
nepochs = 16

"geometry"
# loads design layout
include("layout.jl")
@unpack wwg, hwg, lwg, ld, hclad, hbox, l, w, h, wm, o, base, sources, ports, designs, dx = layout
ϵbox = ϵclad = ϵslab = 2.25f0
ϵcore = 12.25f0

ϵ = sandwich(base, round.(Int, [hbox, hwg, hclad] / dx), [ϵbox, ϵcore, ϵclad])
sz0 = size(ϵ)

"boundaries"
boundaries = [] # unspecified boundaries default to PML

"monitors"
z = [hbox, hbox + hwg]
n = [1, 0, 0]
monitors = vcat([Monitor([x, y, z], n) for (x, y) = ports], monitors_on_box([o..., hbox], [ld, ld, hwg]))


"sources"

@load "$(@__DIR__)/mode.bson" mode
@unpack Ex, Ey, Ez, bounds = mode
f32(x::Real) = Float32(x)
f32(x::AbstractArray) = f32.(x)
bounds = f32(bounds)
bounds /= λ
_dx = mode.dx / λ
Ex, Ey, Ez = [Ex, Ey, Ez] / maximum(maximum.(abs, [Ex, Ey, Ez]))
sy, sz = size.((Ex,), (1, 2))
y = range(F.(bounds[1])..., sy)
z = range(F.(bounds[2])..., sz)

gEx = LinearInterpolation((y, z), F.(Ez))
gEy = LinearInterpolation((y, z), F.(Ex))
gEz = LinearInterpolation((y, z), F.(Ey))
center = [sources[1]..., hbox]
bounds = [0, bounds...]
sources = [
    Source(t -> cos(F(2π) * t), center, bounds;
        Jx=(x, y, z) -> gEx(y, z),
        Jy=(x, y, z) -> gEy(y, z),
        Jz=(x, y, z) -> gEz(y, z),
    ),
]

configs = setup(boundaries, sources, monitors, dx, sz0; F, Courant, T)
@unpack μ, σ, σm, dt, geometry_padding, geometry_splits, field_padding, source_effects, monitor_instances, u0, = configs



function make_geometry(model, μ, σ, σm)

    b = place!(F.(base), model(), round.(Int, o / dx) .+ 1)
    ϵ = sandwich(b, round.(Int, [hbox, hwg, hclad] / dx), [ϵbox, ϵcore, ϵclad])
    ϵ, μ, σ, σm = apply(geometry_padding; ϵ, μ, σ, σm)
    p = apply(geometry_splits; ϵ, μ, σ, σm)
end


# run pre or post optimization simulation and save movie

function metrics(model, T=T;)
    p = make_geometry(model, μ, σ, σm)
    # run simulation
    u = reduce((u, t) -> step!(u, p, t, dx, dt, field_padding, source_effects;), 0:dt:T-1, init=deepcopy(u0))
    reduce(((u, y), t) -> (
            step!(u, p, t, dx, dt, field_padding, source_effects),
            y + dt * power.(monitor_instances, (u,))
        ),
        T-1+dt:dt:T,
        init=(u, zeros(F, length(monitor_instances))))[2]
end

function loss(y)
    dot([0, -1, -1, 1, -1, 1, 1, 1, 1], y) +
    abs(y[2] - y[3])
end

"geometry generator model"
contrast = 10.0f0
nbasis = 4
model = Mask(round.(Int, (ld, ld) ./ dx), nbasis, contrast)#, symmetries=2)
model0 = deepcopy(model)


p0 = make_geometry(model0, μ, σ, σm)
tp = abs(metrics(model, 2)[1])

# volume(p0[1][3])

"Optim functions"
x0, re = destructure(model)
f_ = loss ∘ metrics
f = f_ ∘ re
x = copy(x0)
@show f(x)

error()
# @show x|>re|>metrics

# train surrogate
Random.seed!(1)
runs = []
i = 1
function f!(x)
    @time begin
        model = re(x)
        mp = metrics(model)
        l = loss(mp)
        global i
        # if i == 1
        push!(runs, [x, mp, l])
        # else
        #     runs[i] .= [x0, mp, l]
        # end
        # i +=1
        l
    end
end

iterations = 24
f!(x)
# runs = similar(runs, iterations)

@showtime res = optimize(f!, x0, ParticleSwarm(;
        n_particles=16), Optim.Options(; f_tol=0, iterations, show_every=1, show_trace=true))
xgf = minimizer(res)
x = copy(xgf)

error()
X = Base.stack(getindex.(runs, 1))
Y = Base.stack(getindex.(runs, 2))

nl = leakyrelu
n = size(X, 1)
m, N = size(Y)
nn = Chain(Dense(n, 2n, nl), Dense(2n, 4n, nl), Dense(4n, m))

opt = Adam(0.1)
opt_state = Flux.setup(opt, nn)
n = 100
for i = 1:n
    l, (dldm,) = withgradient(nn -> Flux.mae(Y, nn(X)), nn)
    Flux.update!(opt_state, nn, dldm)
    (i % 25 == 0 || i == n) && println("$i $l")
end

opt = Adam(0.1)
opt_state = Flux.setup(opt, x)
n = 100
@show f(x)
for i = 1:n
    l, (dldm,) = withgradient(x -> mae(
            nn(x)[2:end],
            [tp / 2, tp / 2, -tp, tp, zeros(F, 4)...]
        ), x)
    Flux.update!(opt_state, x, dldm)
    (i % 25 == 0 || i == n) && println("$i $l")
end
@show loss(nn(x))
@show f(x)
xsg = x
x = copy(xsg)

# runsavesim(model0)
# loss(model)
# # volume(p[1][1])

# @showtime sol = accumulate((u, t) -> step(u, p, t, field_padding, source_effects), 0:dt:T, init=u0)

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
    p = make_geometry(model, μ, σ, σm)
    t = 0:dt:T
    sol = similar([u0], length(t))
    @showtime sol = accumulate!((u, t) -> step!(u, p, t, dx, dt, field_padding, source_effects), sol, t, init=u0)

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