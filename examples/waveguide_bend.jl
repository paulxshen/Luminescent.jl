"""
Solves first photonics inverse design problem in Google's Ceviche challenges. 
Design is a tight right angle waveguide bend.
"""

using UnPack, LinearAlgebra, Random, Images, Interpolations, Flux
using Flux: withgradient
using Zygote: ignore, Buffer
using BSON: @save, @load
# using Jello, ArrayPadding

include("../../ArrayPadding.jl/src/ArrayPadding.jl")
include("../../Jello.jl/src/Jello.jl")
using .Jello, .ArrayPadding

p = "../src"
include("$p/startup.jl")
include("$p/del.jl")
include("$p/boundaries.jl")
include("$p/sources.jl")
include("$p/monitors.jl")
include("$p/utils.jl")
include("$p/fdtd.jl")

using CairoMakie
include("$p/plot_recipes.jl")

F = Float32
Random.seed!(1)

# tunable configs
nepochs = 128
T = 10.0f0 # simulation duration in [periods]
opt = Adam(1) # higher learning rate helps
dx = 1.0f0 / 32 # pixel resolution in [wavelengths]
α = 0.2 # grayscale gradient in design mask
lmin = 0.2f0  # minimum feature length in design region [microns]
Courant = 0.5f0 # Courant number
C = 1000 # scale loss

# recording params
frameat = 1 / 4
framerate = 16

# loads Google's Ceviche challenge
λ = 1.28f0 # wavelength in microns
_dx = λ * dx
@unpack base, design_start, design_sz, ports, TE0 = load("waveguide_bend.bson", _dx)
L = size(base) .* dx # domain dimensions [wavelength]
sz0 = size(base)
tspan = (0.0f0, T)
dt = dx * Courant
x, v = TE0
x ./= λ
v /= maximum(v)

# setup FDTD
polarization = :TMz
boundaries = [] # unspecified boundaries default to PML
monitors = [Monitor(ports[1] / λ, [:Ez, :Hy]), Monitor(ports[2] / λ, [:Ez, :Hx,])]
g = linear_interpolation(x, v)
sources = [CenteredSource(t -> cos(F(2π) * t), (x, y) -> g(y), ports[1] / λ, [0, 0.6 / λ]; Jz=1)]
@unpack geometry_padding, field_padding, source_effects, monitor_configs, fields =
    setup(boundaries, sources, monitors, L, dx, polarization; F)

a = ones(F, sz0)
μ = 1a
σ = 0a
σm = 0a
μ, σ, σm = apply(geometry_padding; μ, σ, σm)
u0 = collect(values(fields))
sz = size(first(fields)) # full field size including PML padding

# setup design region to be optimized
model = Jello.Mask(design_sz, lmin / λ / dx; diagonal_symmetry=true) # parameterized binary mask for design region
model0 = deepcopy(model)
ϵ1 = 2.25f0 # oxide
ϵ2 = 12.25f0 # silicon

∇ = Del([dx, dx])
function step(u, p, t)
    ϵ, μ, σ, σm = [[a] for a = p]
    Ez, Hx, Hy, = u
    Jz, = apply(source_effects, t; Jz=0)

    # first update E
    H_ = apply(field_padding, ; Hx=PaddedArray(Hx), Hy=PaddedArray(Hy))
    E = [Ez]
    J = [Jz]
    dEdt = (∇ × H_ - σ * E - J) / ϵ
    E += dEdt * dt
    Ez, = E

    # then update H
    E_ = apply(field_padding, ; Ez=PaddedArray(Ez))
    H = [Hx, Hy]
    dHdt = -(∇ × E_ + σm * H) / μ
    H += dHdt * dt
    Hx, Hy = H

    [Ez, Hx, Hy,]
end

function make_geometry(base, design, start)
    b = place(base, design, start)
    ϵ = ϵ2 * b + ϵ1 * (1 .- b)
    ϵ, = apply(geometry_padding; ϵ)
    p = [ϵ, μ, σ, σm]
end
# pre/post-optimization simulation
function sim(u0, p)
    u = u0
    [
        begin
            u = step(u, p, t)
        end for t = 0:dt:T
    ]
end


function loss(model, α=0.1)
    u = u0
    p = make_geometry(base, model(α), design_start)

    # run first T - 1 periods
    for t = 0:dt:T-1
        u = step(u, p, t)
    end

    # objective to maximize outgoing power at port 2 during last period
    idxs = monitor_configs[2].idxs
    C * reduce(((u, l), t) -> (step(u, p, t), l - (u[1][idxs...])^2), T-1+dt:dt:T, init=(u, 0.0f0))[2]
    # reduce(((u, l), t) -> (step(u, p, t), l + (u[1][idxs...])*u[2][idxs...]), T-1+dt:dt:T, init=(u, 0.0f0))[2]
end
# pre-optimization simulation
p0 = make_geometry(base, model0(0), design_start)
@showtime sol0 = sim(u0, p0)
recordsim(sol0, p0, "bend-pre.gif", title="start of training"; frameat, framerate)

opt_state = Flux.setup(opt, model)
# for i = 1:32

for i = 1:nepochs
    @time begin
        l, (dldm,) = withgradient(m -> loss(m, α), model,)
        Flux.update!(opt_state, model, dldm)
        println("$i $l")
    end
end


if nepochs > 0
    # post-optimization simulation
    p = make_geometry(base, model(0), design_start)
    @showtime sol = sim(u0, p)

    # save movies of simulation
    recordsim(sol, p, "bend-post.gif", title="after some training"; frameat, framerate)
end



plotmonitors(sol0, monitor_configs,)
