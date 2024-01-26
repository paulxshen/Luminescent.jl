using UnPack, LinearAlgebra, Random, StatsBase, Interpolations, Zygote, Optim, Jello
using Zygote: withgradient, Buffer
using Optim: Options, minimizer
using BSON: @save, @load

using GLMakie
dir = pwd()
include("$dir/src/main.jl")
# dir = "FDTDEngine.jl"
# using FDTDEngine

include("$dir/scripts/plot_recipes.jl")
include("$dir/scripts/startup.jl")


F = Float32
Random.seed!(1)

"training params"
name = "silicon_photonics_splitter"
nres = 16
T = 1.0f0 # simulation duration in [periods]
nbasis = 4 # complexity of design region
Courant = 0.5f0 # Courant number
C = 1000 # scale loss
λ = 1.55f0
nepochs = 1

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
    Monitor([[ox, ox + ld], [oy, oy + ld], oz], [0, 0, -1]),
    Monitor([[ox, ox + ld], [oy, oy + ld], oz + hwg], [0, 0, 1]),
    Monitor([[ox, ox + ld], oy, [oz, oz + hwg]], [0, -1, 0]),
    Monitor([[ox, ox + ld], oy + ld, [oz, oz + hwg]], [0, 1, 0]),
    Monitor([ox, [oy, oy + ld], [oz, oz + hwg]], [-1, 0, 0]),
    Monitor([ox + ld, [oy, oy + ld], [oz, oz + hwg]], [1, 0, 0]),
])

"sources"

@load "$(@__DIR__)/mode.bson" mode
@unpack Ex, Ey, Ez, bounds = mode
bounds = f32(bounds)
bounds /= λ
_dx = mode.dx / λ
Ex, Ey, Ez = [Ex, Ey, Ez] / maximum(maximum.(abs, [Ex, Ey, Ez]))
sy, sz = size.((Ex,), (1, 2))
y = range(F.(bounds[1])..., sy)
z = range(F.(bounds[2])..., sz)

f = t -> cos(F(2π) * t)
gEx = LinearInterpolation((y, z), F.(Ez))
gEy = LinearInterpolation((y, z), F.(Ex))
gEz = LinearInterpolation((y, z), F.(Ey))
center = [sources[1]..., hbox]
bounds = [0, bounds...]
sources = [
    Source(f, center, bounds;
        Jx=(x, y, z) -> gEx(y, z),
        Jy=(x, y, z) -> gEy(y, z),
        Jz=(x, y, z) -> gEz(y, z),
    ),
]

configs = setup(boundaries, sources, monitors, dx, sz0; F, Courant, T)
@unpack μ, σ, σm, dt, geometry_padding, geometry_splits, field_padding, source_effects, monitor_instances, fields, step, power = configs



function make_geometry(model, μ, σ, σm)

    b = place(base, model(), round.(Int, o / dx) .+ 1)
    ϵ = sandwich(b, round.(Int, [hbox, hwg, hclad] / dx), [ϵbox, ϵcore, ϵclad])
    ϵ, μ, σ, σm = apply(geometry_padding; ϵ, μ, σ, σm)
    p = apply(geometry_splits; ϵ, μ, σ, σm)
end

p = make_geometry(model, μ, σ, σm)
volume(p[1][1])
# run pre or post optimization simulation and save movie

function savesim(sol, p, name=name)
    E = map(sol) do u
        # sqrt(sum(u) do u u.^2 end)
        u[1] .^ 2
    end
    dir = @__DIR__
    recordsim(E, p[1][1], configs, "$dir/$(name)_nres_$nres.mp4", title="$name"; playback=1, bipolar=false)
end
function loss(model)
    p = make_geometry(model, μ, σ, σm)
    # run simulation
    u = reduce((u, t) -> step(u, p, t, configs; Buffer), 0:dt:T, init=u0)
    reduce(((u, l), t) -> (step(u, p, t, configs; Buffer), begin
            l + dt * (sum(power.(monitor_instances[4:9], (u,))) - 2sum(power.(monitor_instances[2:3], (u,))))
            # l + sum(sum(u))
        end),
        T-1+dt:dt:T, init=(u, 0.0f0))[2]
end
"geometry generator model"
contrast = 10.0f0
nbasis = 4
model = Mask(round.(Int, (ld, ld) ./ dx), nbasis, contrast, symmetries=2)
# @showtime loss(model)
model0 = deepcopy(model)

u0 = collect(values(fields))
p0 = make_geometry(model0, μ, σ, σm)
# volume(p0[3][1])

"Optim functions"
x0, re = destructure(model)
f = loss ∘ re
function g!(storage, x)
    model = re(x)
    g, = gradient(loss, model)
    storage .= realvec(g.a)
end
function fg!(storage, x)
    model = re(x)
    l, (g,) = withgradient(loss, model)
    storage .= realvec(g.a)
    l
end

# runsavesim(model0)
# loss(model)

@showtime sol = accumulate((u, t) -> step(u, p, t, configs), 0:dt:T, init=u0)

