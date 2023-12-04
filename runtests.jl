using CairoMakie
using BSON: @save, @load
using DifferentialEquations, OrdinaryDiffEq, ComponentArrays, UnPack, LinearAlgebra, Random, DiffEqFlux

# using Jello
using Images
include("startup.jl")
include("../Porcupine.jl/src/del.jl")
include("../Jello.jl/src/mask.jl")
F = Float32
include("fdtd.jl")
include("saveimg.jl")

Random.seed!(1)
train = true
# train = false

@load "bend.bson" base design_start design_sz dx
# heatmap(base)
T = 4.0f0
tspan = (0.0f0, T)
microns_dx = dx
dx = 1.0f0 / 16
dt = dx / 2
solver = Euler()

L = size(base) .* microns_dx

λ = 1.28f0
L = L ./ λ
l, _ = L
lc = (57.5f0 / 108) * l
sz = round.(Int, L ./ dx)
design_sz = round.(Int, design_sz .* microns_dx ./ λ ./ dx)
design_start = 1 .+ round.(Int, (design_start .- 1) .* microns_dx ./ λ ./ dx)
base = imresize(base, sz)
model = Mask(design_sz, 0.2f0 / dx)
if train
end
# heatmap(base)
polarization = :TMz
ϵ1 = 2.25f0
ϵ2 = 12.25f0
b = place(base, model(), design_start)
ϵ = ϵ2 * b + ϵ1 * (1 .- b)
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

    u = apply(source_effects, u, u.t)
    E, H, J = group.((u,), [:E, :H, :J])

    H_ = apply(field_padding, (; Hx=u.Hx, Hy=u.Hy);)
    dEzdt, = (∇ × H_ - σ * E - J) / ϵ
    E += [dEzdt * dt]

    E_ = apply(field_padding, (; Ez=E[1]);)
    dHxdt, dHydt = -(∇ × E_ + σm * H) / μ
    dJzdt = zeros(F, size(dEzdt))
    # ComponentArray((; Ez, Hx, Hy, Jz,))
    comp_vec((:Ez, :Hx, :Hy, :Jz, :t), dEzdt, dHxdt, dHydt, dJzdt, 1.0f0,)
end

dudt = Du(p)
u0 = fields
tstops = range(tspan..., length=8 + 1)
callback = train ? nothing : PresetTimeCallback(tstops, saveimg)
prob_neuralode = NeuralODE(dudt, tspan, solver, ; save_idxs, dt, callback, saveat=dt)
function loss(model; withsol=false)
    b = place(base, model(), design_start)
    ϵ = ϵ2 * b + ϵ1 * (1 .- b)
    geometry = (; ϵ, static_geometry...)
    geometry = apply(geometry_padding, geometry)
    # @unpack ϵ, μ, σ, σm = geometry

    # p = comp_vec((:ϵ, :μ, :σ, :σm), Array(ϵ), Array(μ), Array(σ), Array(σm),)
    p = reduce(vcat, vec.(values(geometry)))
    sol = Array(prob_neuralode(u0, p, ;))
    res = make(monitor_configs, sol,)

    n = round(Int, 1 / dt)
    l = res[2].Ez[end-n+1:end] ⋅ res[2].Hx[end-n+1:end]

    withsol ? (l, sol) : l
end


if train
    g = Zygote.gradient(loss, model)[1]
    error()
    opt = Adam(0.1)
    opt_state = Flux.setup(opt, model)

    # fig = Figure()
    # heatmap(fig[1, 1], m(0), axis=(; title="start of training"))
    for i = 1:4
        global sol
        l, (dldm,) = withgradient(loss, model,)
        # (l, sol), (dldm,) = withgradient(loss, model,)
        update!(opt_state, model, dldm)
        println("$i $l")
    end
    # heatmap(fig[1, 2], m(0), axis=(; title="end of training"))
else

    l, sol = loss(model; withsol=true)
end
# display(fig)


using CairoMakie
using CairoMakie: Axis
# t = 0.0f0:dt:T
t = range(0, T, length=size(sol, 2))
res = make(monitor_configs, sol,)
fig = Figure()
for i = 1:2
    ax = Axis(fig[i, 1], title="")
    m = res[i]
    E, H = if i == 1
        m.Ez, m.Hy
    elseif i == 2
        m.Ez, m.Hx
    end
    lines!(ax, t, E)
    lines!(ax, t, H)
    lines!(ax, t, H .* E)
end
display(fig)