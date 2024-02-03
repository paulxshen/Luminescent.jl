using UnPack, LinearAlgebra, Random, StatsBase, Interpolations, Zygote, Jello, Flux
using Flux: mae, Adam
using Zygote: withgradient, Buffer, bufferfrom
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
λ2 = 1.3f0
nepochs = 16
ω2 = λ / λ2

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
monitors = [Monitor([x, y, z], n) for (x, y) = ports]
ox, oy, = o
oz = hbox
append!(monitors, [
    Monitor([ox, [oy, oy + ld], [oz, oz + hwg]], [-1, 0, 0]),
    Monitor([ox + ld, [oy, oy + ld], [oz, oz + hwg]], [1, 0, 0]),
    Monitor([[ox, ox + ld], oy, [oz, oz + hwg]], [0, -1, 0]),
    Monitor([[ox, ox + ld], oy + ld, [oz, oz + hwg]], [0, 1, 0]),
    Monitor([[ox, ox + ld], [oy, oy + ld], oz], [0, 0, -1]),
    Monitor([[ox, ox + ld], [oy, oy + ld], oz + hwg], [0, 0, 1]),
])

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

f1 = t -> cos(F(2π) * t)
f2 = t -> cos(F(2π) * t * λ / λ2)
gEx = LinearInterpolation((y, z), F.(Ez))
gEy = LinearInterpolation((y, z), F.(Ex))
gEz = LinearInterpolation((y, z), F.(Ey))
center = [sources[1]..., hbox]
bounds = [0, bounds...]
sources = [
    Source(t -> f1(t) + f2(t), center, bounds;
        Jx=(x, y, z) -> gEx(y, z),
        Jy=(x, y, z) -> gEy(y, z),
        Jz=(x, y, z) -> gEz(y, z),
    ),
]

configs = setup(boundaries, sources, monitors, dx, sz0; F, Courant, T)
@unpack μ, σ, σm, dt, geometry_padding, geometry_splits, field_padding, source_effects, monitor_instances, fields, power = configs



function make_geometry(model, μ, σ, σm)

    b = place(base, model(), round.(Int, o / dx) .+ 1)
    ϵ = sandwich(b, round.(Int, [hbox, hwg, hclad] / dx), [ϵbox, ϵcore, ϵclad])
    ϵ, μ, σ, σm = apply(geometry_padding; ϵ, μ, σ, σm)
    p = apply(geometry_splits; ϵ, μ, σ, σm)
end


# run pre or post optimization simulation and save movie
function savesim(sol, p, name=name)
    E = map(sol) do u
        # sqrt(sum(u) do u u.^2 end)
        # u[1] .^ 2
        u[1]
    end
    dir = @__DIR__
    recordsim(E, p[1][1], configs, "$dir/$(name)_nres_$nres.mp4", title="$name"; playback=1, bipolar=false)
end
function monitor_powers(model, T=T; bufferfrom=bufferfrom)
    p = make_geometry(model, μ, σ, σm)
    # run simulation
    u = reduce((u, t) -> step!(u, p, t, dx, dt, field_padding, source_effects; bufferfrom), 0:dt:T-1, init=deepcopy(u0))
    reduce(((u, pow), t) -> (step!(u, p, t, dx, dt, field_padding, source_effects; bufferfrom), begin
            y = power.(monitor_instances, (u,))
            pow + dt * vcat(
                cispi(2t) * y, cispi(2t * ω2) * y)
        end),
        T-1+dt:dt:T, init=(u, zeros(ComplexF32, 2length(monitor_instances))))[2]
end

function loss(p)
    p = abs.(p)
    dot([-1, 1, 1, -1, 1, 1, 1, 1], p[2:9]) + dot([1, -1, 1, -1, 1, 1, 1, 1], p[11:18])
end

"geometry generator model"
contrast = 10.0f0
nbasis = 4
model = Mask(round.(Int, (ld, ld) ./ dx), nbasis, contrast)#, symmetries=2)
model0 = deepcopy(model)

u0 = collect(values(fields))
tp = abs(monitor_powers(model, 2)[1])
p0 = make_geometry(model0, μ, σ, σm)
# volume(p0[3][1])

"Optim functions"
x0, re = destructure(model)
f_ = loss ∘ monitor_powers
f = f_ ∘ re
function g!(storage, x)
    model = re(x)
    g, = gradient(f_, model)
    storage .= realvec(g.a)
end
function fg!(storage, x)
    model = re(x)
    l, (g,) = withgradient(f_, model)
    storage .= realvec(g.a)
    l
end
x = copy(x0)
@show f(x)
# @show x|>re|>monitor_powers
# train surrogate
Random.seed!(1)
runs = []
function f!(x)
    @time begin
        model = re(x)
        mp = monitor_powers(model)
        l = loss(mp)
        push!(runs, (x0, mp, l))
        l
    end
end

@showtime res = optimize(f!, x0, ParticleSwarm(;
        n_particles=16), Optim.Options(f_tol=0, iterations=44, show_every=1, show_trace=true))
xgf = minimizer(res)

# error()
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
x = copy(xgf)
opt_state = Flux.setup(opt, x)
n = 100
@show f(x)
for i = 1:n
    l, (dldm,) = withgradient(x -> mae(nn(x),
            [tp / 2, tp / 2, 0, -tp / 2, tp / 2, zeros(F, 4)...,
                tp / 2, 0, tp / 2, -tp / 2, tp / 2, zeros(F, 4)...,]), x)
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


# od = OnceDifferentiable(f, g!, fg!, x0)
# @showtime res = optimize(od, x0, LBFGS(), Optim.Options(f_tol=0, iterations=1, show_every=1, show_trace=true))
# x_adjoint = re(minimizer(res))

model = re(x)
p = make_geometry(model, μ, σ, σm)
t = 0:dt:T
sol = similar([u0], length(t))
@showtime sol = accumulate!((u, t) -> step!(u, p, t, dx, dt, field_padding, source_effects), sol, t, init=u0)

# make movie
Ez = map(sol) do u
    u[3]
end
ϵz = p[1][3]
dir = @__DIR__
° = π / 180
recordsim(Ez, ϵz, configs, "$dir/$(name)_nres_$nres.mp4", title="$name"; elevation=60°, playback=1, bipolar=true)
