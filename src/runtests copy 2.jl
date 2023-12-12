using CairoMakie
using DifferentialEquations, ComponentArrays, UnPack, LinearAlgebra, Random
using NPZ
# using Jello
using Images
include("../Porcupine.jl/src/del.jl")
include("../Jello.jl/src/mask.jl")
include("fdtd.jl")
include("plotstep.jl")

Random.seed!(1)

bg = npzread("bend.npy")
heatmap(bg)
f(x) = 0.1 < x < 0.9
mstart = Tuple(findfirst(f, bg))
msz = 1 .+ Tuple(Tuple(findlast(f, bg)) .- mstart)
dx_ = 1.6 / 40
L = size(bg) .* dx_

λ = 1.28
L = L ./ λ
l, _ = L
lc = (57.5 / 108) * l
dx = 1 / 40
sz = round.(Int, L ./ dx)
msz = round.(Int, msz .* dx_ ./ λ ./ dx)
mstart = round.(Int, mstart .* dx_ ./ λ ./ dx)
bg = imresize(bg, sz)
bg[[a:b for (a, b) = zip(mstart, mstart .+ msz .- 1)]...] .= 0

polarization = :TMz
ϵ1 = 2.25
ϵ2 = 12.25
# l = 32
model = Mask(msz, 0.2 / dx)

∇ = Del([dx, dx])
function f(fields, p, t)
    @unpack ϵ, μ, σ, σm = padded_geometry
    ϵ, μ, σ, σm = [ϵ], [μ], [σ], [σm]
    fields = NamedTuple(fields)
    fields = apply(source_effects, fields, t)
    @unpack E, H, J = fields
    E, H, J = collect.([E, H, J])

    fields = deepcopy(fields)
    fields = apply(field_padding, fields)
    E_ = values(fields.E)
    H_ = values(fields.H)

    dEz, = (∇ × H_ - σ * E - J) / ϵ
    dHx, dHy = -(∇ × E_ + σm * H) / μ
    ComponentArray(
        E=(; z=dEz),
        H=(; x=dHx, y=dHy),
        J=(; z=zeros(size(dEz))),
    )
end

boundaries = []
monitors = [Monitor([0.15l, lc], [:Ez, :Hy]), Monitor([lc, 0.1], [:Ez, :Hy])]
sources = [GaussianBeam(t -> cos(2π * t), 0.15, (0, lc), -1; Ez=1)]
# function loss(model)

b = place(bg, model(), mstart)
ϵ = ϵ2 * b + ϵ1 * (1 .- b)

geometry = (; ϵ)
# sources = [CurrentSource(t->cos(2π*t), (.2, 0.2), (0.2, l / 2); Ez=1)]
@unpack padded_geometry, field_padding, source_effects, monitor_configs, save_idxs, fields =
    setup(geometry, boundaries, sources, monitors, dx, polarization)


u0 = fields
T = 16.0
tspan = (0.0, T)
dt = dx / 2
saveat = 0:dt:T

tstops = range(tspan..., length=8)
@unpack ϵ, μ, σ, σm = padded_geometry
global ϵ, μ, σ, σm
callback = PresetTimeCallback(tstops, plotstep)


prob = ODEProblem(f, u0, tspan)
sol = solve(prob, DP5(); dt, saveat, callback, save_idxs)
t = 0:dt:T
res = make(monitor_configs, sol, t)

Figure()
plot(t, res[1].E.z)
plot(t, res[2].E.z)
display(fig)

# opt = Adam(0.1)
# opt_state = Flux.setup(opt, model)
# data = [[]]

# fig = Figure()
# heatmap(fig[1, 1], m(0), axis=(; title="start of training"))
# for i = 1:100
#     Flux.train!(model, data, opt_state) do m, _
#         l = loss(m)
#         println(l)
#         l
#     end
# end
# heatmap(fig[1, 2], m(0), axis=(; title="end of training"))
# display(fig)