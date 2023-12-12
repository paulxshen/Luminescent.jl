using UnPack, LinearAlgebra, Random, Images
using Zygote
using Zygote: ignore, Buffer
using BSON: @save, @load
using Jello, ArrayPadding

include("startup.jl")
include("del.jl")
include("boundaries.jl")
include("sources.jl")
include("monitors.jl")
include("utils.jl")
include("fdtd.jl")

using CairoMakie
include("plot_recipes.jl")

F = Float32
Random.seed!(1)

# tunable configs
train = true
# train = false
nepochs = 8
T = 10.0f0 # simulation duration in periods
opt = Adam(0.1)
λ = 1.28f0 # wavelength in microns
dx = 1.0f0 / 32 # pixel resolution in wavelengths
lmin = 0.2f0 / λ # minimum feature length in design region in wavelengths
Courant = 0.5f0 # Courant number

# loads Google's Ceviche challenge
_dx = λ * dx
@unpack base, design_start, design_sz, ports = load("bend0", _dx)
L = size(base) .* dx # domain dimensions [wavelength]
sz0 = size(base)
tspan = (0.0f0, T)
dt = dx * Courant

# setup FDTD
polarization = :TMz
boundaries = [] # unspecified boundaries default to PML
monitors = [Monitor(ports[1] / λ, [:Ez, :Hy]), Monitor(ports[2] / λ, [:Ez, :Hx,])]
sources = [GaussianBeam(t -> cos(F(2π) * t), 0.1f0 / λ, ports[1] / λ, 1; Jz=1)]
@unpack geometry_padding, field_padding, source_effects, monitor_configs, fields =
    setup(boundaries, sources, monitors, L, dx, polarization; F)

μ = ones(F, sz0)
σ = zeros(F, sz0)
σm = zeros(F, sz0)
μ, σ, σm = collect.(apply(geometry_padding; μ, σ, σm))
u0 = collect(values(fields))
sz = size(first(fields)) # full field size including PML padding

# setup design region to be optimized
model = Mask(design_sz, lmin / dx) # parameterized binary mask for design region
model0 = deepcopy(model)
ϵ1 = 2.25f0 # oxide
ϵ2 = 12.25f0 # silicon

∇ = Del([dx, dx])
function step(u, p, t)
    ϵ, μ, σ, σm = p |> a -> [[a] for a = a]
    Ez, Hx, Hy, Jz0 = u
    Jz, = apply(source_effects, t; Jz=Jz0)

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

    [Ez, Hx, Hy, Jz0]
end

function make_geometry(model, base, start)
    b = place(base, model(), start)
    ϵ = ϵ2 * b + ϵ1 * (1 .- b)
    ϵ, = apply(geometry_padding; ϵ)
    p = [ϵ, μ, σ, σm]
end

function loss(model)
    u = u0
    p = make_geometry(model, base, design_start)

    for t = 0:dt:T-1
        u = step(u, p, t)
    end
    idxs = monitor_configs[2].idxs
    reduce(((u, l), t) -> (step(u, p, t), l - (u[1][idxs...])^2), T-1+dt:dt:T, init=(u, 0.0f0))[2]
    # reduce(((u, l), t) -> (step(u, p, t), l + (u[1][idxs...])*u[2][idxs...]), T-1+dt:dt:T, init=(u, 0.0f0))[2]
end


if train
    opt_state = Flux.setup(opt, model)
    for i = 1:32
        # for i = 1:nepochs
        @time begin
            l, (dldm,) = withgradient(loss, model,)
            Flux.update!(opt_state, model, dldm)
            println("$i $l")
        end
    end
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
p0 = make_geometry(model0, base, design_start)
p = make_geometry(model, base, design_start)
@showtime sol0 = sim(u0, p0)
@showtime sol = sim(u0, p)

recordsim(sol0, "bend-pre.gif", title="start of training")
recordsim(sol, "bend-post.gif", title="end of training")
plotmonitors(sol, monitor_configs,)
