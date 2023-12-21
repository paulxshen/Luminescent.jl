"""
Solves first photonics inverse design problem in Google's Ceviche challenges. 
Design is a tight right angle waveguide bend.
"""

using UnPack, LinearAlgebra, Random, Images, Interpolations, Flux, Optim
using Flux: withgradient
using Optim: Options, minimizer
using Zygote: ignore, Buffer
using BSON: @save, @load
# using Jello, ArrayPadding

include("../../ArrayPadding.jl/src/utils.jl")
include("../../ArrayPadding.jl/src/pad.jl")
include("../../Jello.jl/src/mask.jl")

p = "../src"
include("$p/startup.jl")
include("$p/del.jl")
include("$p/maxwell.jl")
include("$p/boundaries.jl")
include("$p/sources.jl")
include("$p/monitors.jl")
include("$p/utils.jl")
include("$p/fdtd.jl")

using CairoMakie
include("$p/plot_recipes.jl")

F = Float32
Random.seed!(1)

function make_geometry(design, static_geometry, geometry_padding)
    @unpack start, μ, σ, σm, ϵ1, ϵ2, base = static_geometry
    b = place(base, design, start)
    ϵ = ϵ2 * b + ϵ1 * (1 .- b)
    ϵ, = apply(geometry_padding; ϵ)
    p = [ϵ, μ, σ, σm]
end

# pre/post-optimization simulation
function sim(u0, p, fdtd_configs)
    @unpack dt, T = fdtd_configs
    u = u0
    [u = step_TMz(u, p, t, fdtd_configs) for t = 0:dt:T]
end

function loss(u0, design, static_geometry, fdtd_configs,)
    @unpack dt, T, monitor_configs = fdtd_configs
    u = u0
    p = make_geometry(design, static_geometry, fdtd_configs.geometry_padding)

    # run first T - 1 periods
    for t = 0:dt:T-1
        u = step_TMz(u, p, t, fdtd_configs)
    end

    # objective to maximize outgoing power at port 2 during last period
    idxs = monitor_configs[2].idxs
    reduce(((u, l), t) -> (step_TMz(u, p, t, fdtd_configs), l - (u[1][idxs...])^2), T-1+dt:dt:T, init=(u, 0.0f0))[2]
    # reduce(((u, l), t) -> (step_TMz(u, p, t), l + (u[1][idxs...])*u[2][idxs...]), T-1+dt:dt:T, init=(u, 0.0f0))[2]
end

model = nothing
# schedule = [(8, 1f-2, 100), (16,1f-4, 400), (32,1f-3, 1)]
schedule = [(16, 0.2, 1.0f-2, 50), (16, 0.4, 5.0f-3, 50)]
for (nres, contrast, f_reltol, iterations) in schedule
    # tunable configs
    T = 10.0f0 # simulation duration in [periods]
    opt = Adam(0.2) # higher learning rate helps
    dx = 1.0f0 / nres # pixel resolution in [wavelengths]
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

    x, v = TE0
    x ./= λ
    v /= maximum(v)
    global model
    if isnothing(model)
        model = Mask(design_sz, 3, contrast; diagonal_symmetry=true) # parameterized binary mask for design region
    else
        global model.sz = design_sz
    end
    # setup FDTD
    polarization = :TMz
    boundaries = [] # unspecified boundaries default to PML
    monitors = [Monitor(ports[1] / λ, [:Ez, :Hy]), Monitor(ports[2] / λ, [:Ez, :Hx,])]
    g = linear_interpolation(x, v)
    sources = [CenteredSource(t -> cos(F(2π) * t), (x, y) -> g(y), ports[1] / λ, [0, 0.6 / λ]; Jz=1)]
    fdtd_configs = setup(boundaries, sources, monitors, L, dx, polarization; F, Courant, T)
    @unpack dt, geometry_padding, field_padding, source_effects, monitor_configs, fields = fdtd_configs

    a = ones(F, sz0)
    μ = 1a
    σ = 0a
    σm = 0a
    μ, σ, σm = apply(geometry_padding; μ, σ, σm)

    u0 = collect(values(fields))
    sz = size(first(fields)) # full field size including PML padding

    # setup design region to be optimized
    model0 = deepcopy(model)
    ϵ1 = 2.25f0 # oxide
    ϵ2 = 12.25f0 # silicon
    static_geometry = (; base, start=design_start, ϵ1, ϵ2, μ, σ, σm)

    # pre-optimization simulation
    p0 = make_geometry(model0(), static_geometry, fdtd_configs.geometry_padding)
    @showtime global sol0 = sim(u0, p0, T, fdtd_configs)
    recordsim(sol0, p0, fdtd_configs, "bend-pre_$nres.gif", title="start of training"; frameat, framerate)

    # opt_state = Flux.setup(opt, model)
    # for i = 1:iterations
    #     @time begin
    #         l, (dldm,) = withgradient(m -> loss(u0, m(), static_geometry, fdtd_configs), model,)
    #         Flux.update!(opt_state, model, dldm)
    #         println("$i $l")
    #     end
    # end
    x0, re = destructure(model)
    _loss = m -> C * loss(u0, m(), static_geometry, fdtd_configs)
    f, g!, fg! = optimfuncs(_loss, re)
    od = OnceDifferentiable(f, g!, fg!, x0)
    @showtime res = optimize(od, x0, LBFGS(), Optim.Options(; f_reltol, g_tol=0, iterations, show_every=1, show_trace=true))
    # res = 

    model = re(minimizer(res))

    if iterations > 0
        # post-optimization simulation
        p = make_geometry(model(), static_geometry, fdtd_configs.geometry_padding)
        @showtime global sol = sim(u0, p, T, fdtd_configs)

        # save movies of simulation
        recordsim(sol, p, fdtd_configs, "bend-post_$nres.gif", title="after some training"; frameat, framerate)
    end
end



plotmonitors(sol0, monitor_configs,)
