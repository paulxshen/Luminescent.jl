using CairoMakie
using BSON: @save, @load
using DifferentialEquations, OrdinaryDiffEq, ComponentArrays, UnPack, LinearAlgebra, Random, DiffEqFlux

# using Jello
using Images
include("startup.jl")
include("../Porcupine.jl/src/del.jl")
include("../Jello.jl/src/Base.jl")
F = Float32
include("utils.jl")
include("setup.jl")
include("plotstep.jl")

Random.seed!(1)

train = true
train = false
@load "waveguide_bend.bson" mask design_start design_dims dx
# heatmap(mask)
T = 5.0f0
tspan = (0.0f0, T)
microns_dx = dx
dx = 1.0f0 / 16
dt = dx / 2
saveat = 1 / 16.0f0
solver = Euler()

L = size(mask) .* microns_dx

λ = 1.28f0
L = L ./ λ
l, _ = L
common_left_pad_amount = (57.5f0 / 108) * l
dims = round.(Int, L ./ dx)
design_dims = round.(Int, design_dims .* microns_dx ./ λ ./ dx)
design_start = 1 .+ round.(Int, (design_start .- 1) .* microns_dx ./ λ ./ dx)
mask = imresize(mask, dims)
model = Blob(design_dims, 0.2f0 / dx)
if train
end
# heatmap(mask)
polarization = :TM
ϵ1 = 2.25f0
ϵ2 = 12.25f0
b = place(mask, model(), design_start)
ϵ = ϵ2 * b + ϵ1 * (1 .- b)
# extrema(ϵ)

boundaries = []
monitors = [Monitor([0.0f0l, common_left_pad_amount], [:Ez, :Hx, :Hy]), Monitor([common_left_pad_amount, 0.0f0], [:Ez, :Hx, :Hy])]
sources = [GaussianBeam(t -> cos(F(2π) * t), 0.05f0, (0.0f0, common_left_pad_amount), -1; Ez=1)]
@unpack geometry_padding, field_padding, source_instances, monitor_configs, save_idxs, fields =
    setup(boundaries, sources, monitors, L, dx, polarization; F)

static_geometry = (; μ=ones(Int, dims), σ=zeros(Int, dims), m=zeros(Int, dims))
geometry = (; ϵ, static_geometry...)
geometry = apply(geometry_padding, geometry)
@unpack ϵ, μ, σ, m = geometry
_dims = size(ϵ)

p = comp_vec((:ϵ, :μ, :σ, :m), Array.((ϵ, μ, σ, m))...,)

∇ = Del([dx, dx])
function dudt(u, p, t)
    ϵ, μ, σ, m = eachslice(reshape(p, _dims..., 4), dims=3)
    # @unpack ϵ, μ, σ, m=m.p
    ϵ, μ, σ, m = [ϵ], [μ], [σ], [m]

    u = apply(source_instances, u, u.t)
    E, H, J = group.((u,), [:E, :H, :J])

    H_ = apply(field_padding, (; Hx=u.Hx, Hy=u.Hy);)
    dEzdt, = (∇ × H_) / ϵ
    # dEzdt, = (∇ × H_ - σ * E - J) / ϵ
    E += [dEzdt * dt]

    E_ = apply(field_padding, (; Ez=E[1]);)
    dHxdt, dHydt = -(∇ × E_ + m * H) / μ
    dJzdt = zeros(Int, size(dEzdt))
    # ComponentArray((; Ez, Hx, Hy, Jz,))
    comp_vec((:Ez, :Hx, :Hy, :Jz, :t), dEzdt, dHxdt, dHydt, dJzdt, 1.0f0,)
end

fields = ComponentArray(merge(u0, (; t=F(0))))
u0 = fields
tstops = range(tspan..., length=8 + 1)
callback = train ? nothing : PresetTimeCallback(tstops, plotstep)
prob_neuralode = ODEProblem(dudt, u0, tspan, p,)
function loss(model; withsol=false)
    b = place(mask, model(), design_start)
    ϵ = ϵ2 * b + ϵ1 * (1 .- b)
    geometry = (; ϵ, static_geometry...)
    geometry = apply(geometry_padding, geometry)
    # @unpack ϵ, μ, σ, m = geometry

    # p = comp_vec((:ϵ, :μ, :σ, :m), Array(ϵ), Array(μ), Array(σ), Array(m),)
    p = reduce(vcat, vec.(values(geometry)))
    # global sol = Array(solve(prob_neuralode, solver; dt, p,
    sol = solve(
        prob_neuralode,
        solver;
        dt,
        p,
        callback,
        saveat,
        sensealg=QuadratureAdjoint(autojacvec=ZygoteVJP(), abstol=1e-2,
            reltol=1e-1),
    )

    t = (round(Int, (T - 1) / saveat):round(Int, T / saveat)) .+ 1
    res = field(Array(sol), monitor_configs[2], t)
    # t = T-1:dt:T
    # res = field(sol, monitor_configs, t)

    l = 1mean(res.Ez .* res.Hx)

    withsol ? (l, sol) : l
end


if train
    # g = Zygote.gradient(loss, model)[1]
    # error()
    opt = Adam(0.1)
    # opt = AdaGrad()
    opt_state = Flux.setup(opt, model)

    # fig = Figure()
    # heatmap(fig[1, 1], m(0), axis=(; title="start of training"))
    for i = 1:2
        global sol, l, dldm
        l, (dldm,) = withgradient(loss, model,)
        # (l, sol), (dldm,) = withgradient(loss, model,)
        Flux.update!(opt_state, model, dldm)
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
fig = Figure()
for i = 1:2
    ax = Axis(fig[i, 1], title="")
    m = field(sol, monitor_configs[i], :)
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