using CairoMakie
using DifferentialEquations, ComponentArrays, UnPack, LinearAlgebra, Random
using NPZ
using Jello
using Images
include("../Porcupine.jl/src/del.jl")
include("fdtd.jl")
Random.seed!(1)

bg = npzread("bend.npy")
heatmap(bg)
f(x) = 0.1 < x < 0.9
mstart = findfirst(f, bg)
msz = 1 .+ collect(findlast(f, bg) - mstart)
dx = 1.6 / 40

λ = 1.28
sz = round.(Int, size(bg) .* dx / λ)
bg = resize(bg, dz)

dx = 1 / 16
polarization = :TMz
ϵ1 = 2.25
ϵ2 = 12.25
# l = 32
# geometry = (; ϵ=ones(l, l))
mask = Mask(msz, 0.2 / dx)
mask = place(bg, mask, mstart)
ϵ = ϵ2 * mask + ϵ1 * (1 .- mask)
geometry = (; ϵ)
boundaries = [Periodic(2), PEC(1)]
sources = [PlaneWave(cos, (; E=(; z=1)), -1)]
@unpack padded_geometry, boundary_effects, source_effects, fields = setup(geometry, boundaries, sources, dx, polarization)


∇ = Del([dx, dx])
function f(fields, p, t)
    @unpack ϵ, μ, σ, σm = padded_geometry
    ϵ, μ, σ, σm = [ϵ], [μ], [σ], [σm]
    fields = NamedTuple(fields)
    fields = apply(source_effects, fields, t)
    @unpack E, H, J = fields

    fields = apply(Array, deepcopy(fields))
    fields = apply(boundary_effects, fields)
    E_ = values(fields.E)
    H_ = values(fields.H)

    dEz, = (∇ × H_ - σ * E - J) / ϵ
    dHx, dHy = -(∇ × E_ + σm * H) / μ
    # dBx, dBy = -∇ × E
    # dHx, dHy = dBx ./ μ, dBy ./ μ

    dE = (; z=dEz)
    dH = (; x=dHx, y=dHy)
    ComponentArray(E=dE, H=dH)
end
u0 = ComponentArray(fields)
T = 1.0
tspan = (0.0, T)
prob = ODEProblem(f, u0, tspan)
saveat = range(0, T, 16)
sol = solve(prob, Euler(); dt=dx / 4, saveat)

a = sol(T)
fig = heatmap(a.E.z)

for t = saveat
    a = sol(t)
    fig = heatmap(a.E.z)
    save("$t.png", fig)
end

@unpack ϵ, μ, σ, σm = padded_geometry
heatmap(σ)