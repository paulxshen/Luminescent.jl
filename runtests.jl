using CairoMakie
using DifferentialEquations, ComponentArrays, UnPack, LinearAlgebra, Random, DiffEqFlux

# using Jello
using Images
include("startup.jl")
include("../Porcupine.jl/src/del.jl")
include("../Jello.jl/src/mask.jl")
F = Float32
include("fdtd.jl")
include("saveimg.jl")
include("invrs.jl")

Random.seed!(1)

train = true

T = 14.0f0
tspan = (0.0f0, T)
dx = 1.0f0 / 40
dt = dx / 2

base, design_start, design_sz, microns_dx = invrs_load("bend")
L = size(base) .* microns_dx

λ = 1.28f0
L = L ./ λ
l, _ = L
lc = (57.5f0 / 108) * l
sz = round.(Int, L ./ dx)
design_sz = round.(Int, design_sz .* microns_dx ./ λ ./ dx)
design_start = round.(Int, design_start .* microns_dx ./ λ ./ dx)
base = imresize(base, sz)
base[[a:b for (a, b) = zip(design_start, design_start .+ design_sz .- 1)]...] .= 0
if train
    model = Mask(design_sz, 0.2f0 / dx)
end

polarization = :TMz
ϵ1 = 2.25f0
ϵ2 = 12.25f0
ϵ = ϵ2 * base + ϵ1 * (1 .- base)
geometry = (; ϵ)
boundaries = []
monitors = [Monitor([0.15f0l, lc], [:Ez, :Hx, :Hy]), Monitor([lc, 0.1f0], [:Ez, :Hx, :Hy])]
sources = [GaussianBeam(t -> cos(F(2π) * t), 0.05f0, (0.0f0, lc), -1; Ez=1)]
@unpack geometry_padding, field_padding, source_effects, monitor_configs, save_idxs, fields =
    setup(boundaries, sources, monitors, L, dx, polarization; F)

static_geometry = (; μ=ones(F, sz), σ=zeros(F, sz), σm=zeros(F, sz))
geometry = (; ϵ, static_geometry...)
geometry = apply(geometry_padding, geometry)
@unpack ϵ, μ, σ, σm = geometry
_sz = size(ϵ)

p = comp_vec((:ϵ, :μ, :σ, :σm), Array.((ϵ, μ, σ, σm))...,)

struct Du
    p
end
Flux.@functor Du
∇ = Del([dx, dx])
function (m::Du)(u)
    ϵ, μ, σ, σm = eachslice(reshape(m.p, _sz..., 4), dims=3)
    # @unpack ϵ, μ, σ, σm=m.p
    ϵ, μ, σ, σm = [ϵ], [μ], [σ], [σm]

    u = NamedTuple(u)
    u = apply(source_effects, u, u.t)
    E, H, J = group.((u,), [:E, :H, :J])
    u = apply(field_padding, u)
    E_, H_ = group.((u,), [:E, :H,])

    Ez, = (∇ × H_ - σ * E - J) / ϵ
    Hx, Hy = -(∇ × E_ + σm * H) / μ
    Jz = zeros(F, size(Ez))
    # ComponentArray((; Ez, Hx, Hy, Jz,))
    comp_vec((:Ez, :Hx, :Hy, :Jz, :t), Ez, Hx, Hy, Jz, 1.0f0,)
end

dudt = Du(p)
u0 = fields
tstops = range(tspan..., length=8 + 1)
callback = train ? nothing : PresetTimeCallback(tstops, saveimg)
prob_neuralode = NeuralODE(dudt, tspan, Tsit5(), ; save_idxs, dt, callback, saveat=dt)
function loss(model)
    b = place(base, model(), design_start)
    ϵ = ϵ2 * b + ϵ1 * (1 .- b)
    geometry = (; ϵ, static_geometry...)
    geometry = apply(geometry_padding, geometry)
    @unpack ϵ, μ, σ, σm = geometry

    p = comp_vec((:ϵ, :μ, :σ, :σm), Array(ϵ), Array(μ), Array(σ), Array(σm),)
    global sol = Array(prob_neuralode(u0, p, ;))
    res = make(monitor_configs, sol,)

    n = round(Int, 1 / dt)
    res[2].Ez[end-n+1:end] ⋅ res[2].Hx[end-n+1:end]
end


if train
    opt = Adam(0.1)
    opt_state = Flux.setup(opt, model)
    data = [[]]

    # fig = Figure()
    # heatmap(fig[1, 1], m(0), axis=(; title="start of training"))
    for i = 1:1
        l, (dldm,) = withgradient(loss, model,)
        update!(opt_state, model, dldm)
        println("$i $l")
    end
    # heatmap(fig[1, 2], m(0), axis=(; title="end of training"))
else

    loss(model)
end
# display(fig)


using CairoMakie
using CairoMakie: Axis
# t = 0.0f0:dt:T
t = range(0, T, length=size(sol, 2))
res = make(monitor_configs, sol,)
fig = Figure()
ax = Axis(fig[1, 1], title="")
lines!(ax, res[1].Ez)
lines!(ax, res[2].Hx)
display(fig)