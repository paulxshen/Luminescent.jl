using CairoMakie
using DifferentialEquations, ComponentArrays, UnPack, LinearAlgebra, Random, DiffEqFlux

# using Jello
using Images
include("startup.jl")
include("../Porcupine.jl/src/del.jl")
include("../Jello.jl/src/mask.jl")
include("fdtd.jl")
include("saveimg.jl")
include("invrs.jl")

Random.seed!(1)

base, design_start, design_sz, microns_dx = invrs_load("bend")
L = size(base) .* microns_dx

λ = 1.28
L = L ./ λ
l, _ = L
lc = (57.5 / 108) * l
dx = 1 / 8
sz = round.(Int, L ./ dx)
design_sz = round.(Int, design_sz .* microns_dx ./ λ ./ dx)
design_start = round.(Int, design_start .* microns_dx ./ λ ./ dx)
base = imresize(base, sz)
base[[a:b for (a, b) = zip(design_start, design_start .+ design_sz .- 1)]...] .= 0

polarization = :TMz
ϵ1 = 2.25
ϵ2 = 12.25

model = Mask(design_sz, 0.2 / dx)
T = 16.0
tspan = (0.0, T)
dt = dx / 2
saveat = 0:dt:T
ϵ = ϵ2 * base + ϵ1 * (1 .- base)
geometry = (; ϵ)
boundaries = []
monitors = [Monitor([0.15l, lc], [:Ez, :Hx, :Hy]), Monitor([lc, 0.1], [:Ez, :Hx, :Hy])]
sources = [GaussianBeam(t -> cos(2π * t), 0.15, (0, lc), -1; Ez=1)]
@unpack geometry_effects, boundary_effects, source_effects, monitor_configs, save_idxs, fields =
    setup(boundaries, sources, monitors, L, dx, polarization)

∇ = Del([dx, dx])
function dudt(fields, p, t)
    ϵ, μ, σ, σm = [[p.ϵ], [p.μ], [p.σ], [p.σm]]
    # ϵ, μ, σ, σm = [:ϵ, :μ, :σ, :σm] |> k -> [p[k]]
    fields = NamedTuple(fields)
    fields = apply(source_effects, fields, t)
    E, H, J = group.((fields,), [:E, :H, :J])

    fields = apply(boundary_effects, deepcopy(fields))
    global E_, H_ = group.((fields,), [:E, :H,])

    Ez, = (∇ × H_ - σ * E - J) / ϵ
    Hx, Hy = -(∇ × E_ + σm * H) / μ
    Jz = zeros(size(Ez))
    ComponentArray(; Ez, Hx, Hy, Jz,)
end

# prob_neuralode = NeuralODE(dudt, tspan, Tsit5(); saveat)
sz = size(ϵ)
geometry = ComponentArray(; ϵ, μ=ones(sz), σ=zeros(sz), σm=zeros(sz))
prob = ODEProblem(dudt, fields, tspan, geometry;)
function loss(model)
    b = place(base, model(), design_start)
    ϵ = ϵ2 * b + ϵ1 * (1 .- b)
    geometry = (; ϵ, μ=ones(sz), σ=zeros(sz), σm=zeros(sz))
    # geometry=merge(geometry,(;ϵ))
    geometry = ComponentArray(apply(geometry_effects, geometry))

    _prob = remake(prob; p=geometry)
    sol = solve(prob_, DP5(); save_idxs, saveat, sensealg=SensitivityADPassThrough(), dense=false)

    t = 0:dt:T
    res = make(monitor_configs, sol, t)

    n = round(Int, 1 / dt)
    res[2].Ez[end-n+1:end] ⋅ res[2].Hx[end-n+1:end]
end
tstops = range(tspan..., length=8)
callback = PresetTimeCallback(tstops, saveimg)

loss(model)

opt = Adam(0.1)
opt_state = Flux.setup(opt, model)
data = [[]]

# fig = Figure()
# heatmap(fig[1, 1], m(0), axis=(; title="start of training"))
for i = 1:1
    Flux.train!(model, data, opt_state) do m, _
        l = loss(m)
        println(l)
        l
    end
end
# heatmap(fig[1, 2], m(0), axis=(; title="end of training"))
# display(fig)