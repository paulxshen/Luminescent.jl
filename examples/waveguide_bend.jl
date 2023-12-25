"""
Solves first photonics inverse design problem in Google's Ceviche challenges. 
Design is a tight right angle waveguide bend.
"""

using UnPack, LinearAlgebra, Random, Images, Interpolations, Flux, Optim, DataStructures
using Flux: withgradient
using Optim: Options, minimizer
using Zygote: ignore, Buffer
using BSON: @save, @load
# using Jello, ArrayPadding

include("../../ArrayPadding.jl/src/utils.jl")
include("../../ArrayPadding.jl/src/pad.jl")
include("../../Jello.jl/src/mask.jl")

dir = "../src"
include("$dir/startup.jl")
include("$dir/del.jl")
include("$dir/maxwell.jl")
include("$dir/boundaries.jl")
include("$dir/sources.jl")
include("$dir/monitors.jl")
include("$dir/utils.jl")
include("$dir/fdtd.jl")

using CairoMakie
include("$dir/plot_recipes.jl")

F = Float32
Random.seed!(1)
name = "waveguide_bend"

function make_geometry(design, static_geometry, geometry_padding)
    @unpack start, μ, σ, σm, ϵ1, ϵ2, base = static_geometry
    b = place(base, design, start)
    ϵ = ϵ2 * b + ϵ1 * (1 .- b)
    ϵ = apply(geometry_padding[:ϵ], ϵ)
    p = [ϵ, μ, σ, σm]
end

# pre/post-optimization simulation
function sim(u0, p, fdtd_configs)
    @unpack dt, T = fdtd_configs
    u = u0
    [u = stepTMz(u, p, t, fdtd_configs) for t = 0:dt:T]
end

function loss(u0, design, static_geometry, fdtd_configs,)
    @unpack dx, dt, T, monitor_info = fdtd_configs
    u = u0
    p = make_geometry(design, static_geometry, fdtd_configs.geometry_padding)

    # run first T - 1 periods
    for t = 0:dt:T-1
        u = stepTMz(u, p, t, fdtd_configs)
    end

    # objective to maximize outgoing power at port 2 during last period
    # reduce(((u, l), t) -> (stepTMz(u, p, t, fdtd_configs), l - (u[1][idxs...])^2), T-1+dt:dt:T, init=(u, 0.0f0))[2]
    reduce(((u, l), t) -> (stepTMz(u, p, t, fdtd_configs), begin
            for i = 3:7
                Ez, Hx, Hy = getindex.(u, Ref.(monitor_info[i].idxs)...)
                l +=
                    if i == 2
                        sum(Ez .* Hx + 0 * Hy) # must include all fields otherwise Zygote fails
                    elseif i == 3
                        sum(Ez .* Hy + 0 * Hx)
                    elseif i == 4
                        -sum(Ez .* Hy + 0 * Hx)
                    elseif i == 5
                        sum(Ez .* Hx + 0 * Hy)
                    elseif i == 6 || i == 7
                        -sum(Ez .* Hx + 0 * Hy)
                    end
            end
            l
        end),
        T-1+dt:dt:T, init=(u, 0.0f0))[2] * dx
end

model = nothing
schedule = [(16, 0.9, 5.0f-3, 50)]
# schedule = [(16, 0.5, 1.0f-2, 50), (16, 0.8, 5.0f-3, 50)]
for (nres, contrast, f_reltol, iterations) in schedule
    # tunable configs
    T = 10.0f0 # simulation duration in [periods]
    dx = 1.0f0 / nres # pixel resolution in [wavelengths]
    nbasis = 4 # complexity of design region
    doflux = false
    η = 0.2
    Courant = 0.5f0 # Courant number
    C = 1000 # scale loss

    # recording params
    frameat = 1 / 4
    framerate = 16

    # loads Google's Ceviche challenge
    λ = 1.28f0 # wavelength in microns
    _dx = λ * dx
    @unpack base, design_start, design_sz, wwg, ports, TE0, wwg, lwg, ld, lp = load("examples", "$name.bson", _dx)
    L = size(base) .* dx # domain dimensions [wavelength]
    sz0 = size(base)
    tspan = (0.0f0, T)

    x, v = TE0
    x ./= λ
    v /= maximum(v)

    global model
    if isnothing(model)
        model = Mask(design_sz, nbasis, contrast; diagonal_symmetry=true) # parameterized binary mask for design region
    else
        model = Mask(model; dims=design_sz)
    end

    # setup FDTD
    polarization = :TMz
    boundaries = [] # unspecified boundaries default to PML
    r = 1.6wwg / 2
    monitors = Monitor.([
        [lwg / 2, [lp - r, lp + r]],
        [[lp - r, lp + r], lwg / 2],
        [lwg, [lwg, lwg + ld]],
        [lwg + ld, [lwg, lwg + ld]],
        [[lwg, lwg + ld], lwg + ld,],
        [[lwg, lp - r], lwg,],
        [[lp + r, lwg + ld], lwg,],
        # [lwg + ld, lwg,lwg+ld,],
    ] ./ λ)
    g = linear_interpolation(x, v)
    sources = [CenteredSource(t -> cos(F(2π) * t), (x, y) -> g(y), ports[1] / λ, [0, 0.6 / λ]; Jz=1)]
    global fdtd_configs = setup(boundaries, sources, monitors, L, dx, polarization; F, Courant, T)
    @unpack dt, geometry_padding, field_padding, source_effects, monitor_info, fields = fdtd_configs

    μ = apply(geometry_padding[:μ], ones(F, sz0))
    σ = apply(geometry_padding[:σ], zeros(F, sz0))
    σm = apply(geometry_padding[:σm], zeros(F, sz0))

    u0 = collect(values(fields))
    dims = size(first(fields)) # full field size including PML padding

    # setup design region to be optimized
    model0 = deepcopy(model)
    ϵ1 = 2.25f0 # oxide
    ϵ2 = 12.25f0 # silicon
    static_geometry = (; base, start=design_start, ϵ1, ϵ2, μ, σ, σm)

    # pre-optimization simulation
    p0 = make_geometry(model0(), static_geometry, fdtd_configs.geometry_padding)
    @showtime global sol0 = sim(u0, p0, fdtd_configs)
    recordsim(sol0, p0, fdtd_configs, "$name-pre_$nres.gif", title="start of training"; frameat, framerate)
    if doflux
        "Flux.jl train"
        opt = Adam(η)
        opt_state = Flux.setup(opt, model)
        for i = 1:iterations
            @time begin
                l, (dldm,) = withgradient(model -> C * loss(u0, model(), static_geometry, fdtd_configs), model)
                # l, dldm = withgradient(() -> loss(u0, model(), static_geometry, fdtd_configs), params(model),)
                Flux.update!(opt_state, model, dldm)
                println("$i $l")
            end
        end
    else
        "Optim.jl train"
        x0, re = destructure(model)
        _loss = m -> C * loss(u0, m(), static_geometry, fdtd_configs)
        f, g!, fg! = optimfuncs(_loss, re)
        od = OnceDifferentiable(f, g!, fg!, x0)
        @showtime res = optimize(od, x0, LBFGS(), Optim.Options(; f_reltol, g_tol=0, iterations, show_every=1, show_trace=true))
        model = re(minimizer(res))
    end
    if iterations > 0
        # post-optimization simulation
        p = make_geometry(model(), static_geometry, fdtd_configs.geometry_padding)
        @showtime global sol = sim(u0, p, fdtd_configs)

        # save movies of simulation
        recordsim(sol, p, fdtd_configs, "$name-post_$nres.gif", title="after some training"; frameat, framerate)
    end
end



# plotmonitors(sol0, monitor_info,)
