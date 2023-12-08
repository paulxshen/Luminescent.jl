using CairoMakie
using BSON: @save, @load
using ComponentArrays, UnPack, LinearAlgebra, Random, LazyStack

# using Jello
using Images
include("startup.jl")
include("../Porcupine.jl/src/del.jl")
include("../Jello.jl/src/mask.jl")
F = Float32
include("utils.jl")
include("fdtd.jl")
include("saveimg.jl")

Random.seed!(1)

train = true
# train = false
nepochs = 8
T = 8.0f0
opt = Adam(0.01)

λ = 1.28f0
dx = 1.0f0 / 32
_dx = λ * dx
@unpack base, design_start, design_sz, l = load("bend0", _dx)
L = size(base) .* dx
l = l ./ λ
sz = size(base)
# heatmap(base)
tspan = (0.0f0, T)
dt = dx / 2
saveat = 1 / 16.0f0
callbackat = 1.0f0

model = Mask(design_sz, 0.2f0 / dx)
if train
end
polarization = :TMz
ϵ1 = 2.25f0
ϵ2 = 12.25f0
b = place(base, model(), design_start)
ϵ = ϵ2 * b + ϵ1 * (1 .- b)
# extrema(ϵ)

boundaries = []
monitors = [Monitor([0.0f0, l], [:Ez, :Hx, :Hy]), Monitor([l, 0.0f0], [:Ez, :Hx, :Hy])]
sources = [GaussianBeam(t -> cos(F(2π) * t), 0.05f0, (0.0f0, l), -1; Ez=1)]
@unpack geometry_padding, field_padding, source_effects, monitor_configs, save_idxs, fields =
    setup(boundaries, sources, monitors, L, dx, polarization; F)

static_geometry = (; μ=ones(F, sz), σ=zeros(F, sz), σm=zeros(F, sz))

# Base.eachslice(a) = eachslice(a, dims=ndims(a))
∇ = Del([dx, dx])
function step(u, p, t)
    ϵ, μ, σ, σm = eachslice(p, dims=ndims(p)) |> a -> [[a] for a = a]

    Ez, Hx, Hy, Jz = apply(source_effects, (:Ez, :Hx, :Hy, :Jz), eachslice(u, dims=ndims(u)), t)

    H_ = apply(field_padding, (:Hx, :Hy), PaddedArray.((Hx, Hy)))
    E = [Ez]
    J = [Jz]
    dEdt = (∇ × H_ - σ * E - J) / ϵ
    E += dEdt * dt
    Ez, = E

    E_ = apply(field_padding, (:Ez,), PaddedArray.((Ez,)))
    H = [Hx, Hy]
    dHdt = -(∇ × E_ + σm * H) / μ
    H += dHdt * dt
    Hx, Hy = H
    lazystack([Ez, Hx, Hy, Jz])
end

callback = train ? nothing : saveimg
dsaveat = round(Int, saveat / dt)
dcallbackat = round(Int, callbackat / dt)
function loss(model; withsol=false, callback=nothing)
    b = place(base, model(), design_start)
    ϵ = ϵ2 * b + ϵ1 * (1 .- b)
    geometry = (ϵ, values(static_geometry)...)
    geometry = apply(geometry_padding, (:ϵ, :μ, :σ, :σm), geometry)

    u = lazystack(values(fields))
    sol = [u]
    p = lazystack(values(geometry))
    for (i, t) = enumerate(0:dt:T)
        u = step(u, p, t)
        if mod1(i, dsaveat) == 1
            sol = vcat(sol, [u])
        end
        if mod1(i, dcallbackat) == 1 && !isnothing(callback)
            callback(u, p, t)
        end
    end

    # # sol = reduce(hcat, reduce.(vcat, values.(sol)))
    t = (round(Int, (T - 1) / saveat):round(Int, T / saveat)) .+ 1
    sol = lazystack(Array.(sol))
    Ez, Hx = get(sol, monitor_configs[2], t)

    l = 1mean(Ez .* Hx)
    # maximum(sol[end].Ez)
    withsol ? (l, sol) : l
end


if train
    opt_state = Flux.setup(opt, model)

    # fig = Figure()
    # heatmap(fig[1, 1], m(0), axis=(; title="start of training"))
    for i = 1:2nepochs
        @time begin

            global dldm
            l, (dldm,) = withgradient(loss, model,)
            # (l, sol), (dldm,) = withgradient(loss, model,)
            Flux.update!(opt_state, model, dldm)
            println("$i $l")
        end
    end
    # heatmap(fig[1, 2], m(0), axis=(; title="end of training"))
end
train = false
@showtime l, sol = loss(model; withsol=true, callback)


using CairoMakie
using CairoMakie: Axis
# t = 0.0f0:dt:T
t = range(0, T, length=size(sol)[end])
fig = Figure()
for i = 1:2
    ax = Axis(fig[i, 1], title="")
    E, H = get(sol, monitor_configs[i],)
    #  = if i == 1
    #     m.Ez, m.Hy
    # elseif i == 2
    #     m.Ez, m.Hx
    # end
    lines!(ax, t, E)
    lines!(ax, t, H)
    lines!(ax, t, H .* E)
end
display(fig)