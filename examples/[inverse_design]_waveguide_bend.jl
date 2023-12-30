"""
Solves first photonics inverse design problem in Google's Ceviche challenges. 
Design is a tight right angle waveguide bend.
"""

######################## 
# tunable configs #
###################
name = "waveguide_bend"
data_dir = "examples/utils"

"training params"
T = 10.0f0 # simulation duration in [periods]
nbasis = 4 # complexity of design region
trainer = :Optim # :Flux or :Optim
η = 0.2 # training rate (only applies if Flux is trainer)
Courant = 0.5f0 # Courant number
C = 1000 # scale loss

"""
training schedule
    - nres: resolution per wavelength
    - contrast: edge sharpness in design region
    - f_reltol: objective relative tolerance 
    - iterations 
"""
schedule = [(16, 0.9, 5.0f-3, 50),] # 

"recording params"
frameat = 1 / 4 # captures frame every _ periods
framerate = 16 # playback speed
########################

using UnPack, LinearAlgebra, Random, StatsBase, ImageTransformations, Interpolations, Flux, Optim, Jello
using Flux: withgradient
using Optim: Options, minimizer

include("../src/fdtd_prerelease.jl")
using .fdtd_prerelease

include("utils/invrs_load.jl")
include("utils/plot_recipes.jl")

F = Float32
Random.seed!(1)


model = nothing
for (nres, contrast, f_reltol, iterations) in schedule
    # loads Google's Ceviche challenge at specified resolution 
    λ = 1.28f0 # wavelength in [microns]
    dx = 1.0f0 / nres # pixel resolution in [wavelengths]
    _dx = λ * dx
    @unpack base, design_start, design_sz, wwg, ports, TE0, wwg, lwg, ld, lp = load(data_dir, "$name.bson", _dx)
    L = size(base) .* dx # domain dimensions [wavelength]
    tspan = (0.0f0, T)


    global model
    if isnothing(model)
        model = Mask(design_sz, nbasis, contrast; diagonal_symmetry=true) # parameterized binary mask for design region
    else
        model = Mask(model; dims=design_sz)
    end

    # setup monitor surfaces 
    """
    wwg: waveguide width [wavelengths]
    lwg: waveguide length 
    ld: design region length 
    lp: edge to port center distance 
    """
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

    # loads modal source profile
    x, v = TE0
    x ./= λ
    v /= maximum(v)
    g = linear_interpolation(x, v)
    sources = [CenteredSource(t -> cos(F(2π) * t), (x, y) -> g(y), [lwg / 2, lp] / λ, [0, 0.6 / λ]; Jz=1)]


    # setup FDTD
    polarization = :TMz # Ez, Hx, Hy
    boundaries = [] # unspecified boundaries default to PML
    global fdtd_configs = setup(boundaries, sources, monitors, L, dx, polarization; F, Courant, T)
    @unpack Ie, Ih, dx, esz0, hsz0, esz, hsz, dt, geometry_padding, field_padding, source_effects, monitor_instances, fields = fdtd_configs

    μ = apply(geometry_padding[:μ], ones(F, hsz0))
    σ = apply(geometry_padding[:σ], zeros(F, esz0))
    σm = apply(geometry_padding[:σm], zeros(F, hsz0))

    # esz, hsz: full field size including PML padding

    # setup design region to be optimized
    model0 = deepcopy(model)
    ϵ1 = 2.25f0 # oxide
    ϵ2 = 12.25f0 # silicon
    static_geometry = (; base, start=design_start, ϵ1, ϵ2, μ, σ, σm)

    "places design region onto base geometry"
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

    """
    Objective is to maximize outgoing power at port 2 and minimize power leaked elsewhere  during last period
    The latter auxiliary objectives reduce local minima
    """
    function loss(u0, design, static_geometry, fdtd_configs,)
        @unpack dx, dt, T, monitor_instances = fdtd_configs
        u = u0
        p = make_geometry(design, static_geometry, fdtd_configs.geometry_padding)

        # run first T - 1 periods
        for t = 0:dt:T-1
            u = stepTMz(u, p, t, fdtd_configs)
        end
        reduce(((u, l), t) -> (stepTMz(u, p, t, fdtd_configs), begin
                for i = 3:7
                    Ez, Hx, Hy = getindex.(u, Ref.(monitor_instances[i].idxs)...)
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


    # pre-optimization simulation
    u0 = collect(values(fields))
    p0 = make_geometry(model0(), static_geometry, fdtd_configs.geometry_padding)
    @showtime global sol0 = sim(u0, p0, fdtd_configs)
    recordsim(sol0, p0, fdtd_configs, "$name-pre_$nres.gif", title="start of training"; frameat, framerate)

    if trainer == :Flux
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

    # post-optimization simulation
    if iterations > 0
        p = make_geometry(model(), static_geometry, fdtd_configs.geometry_padding)
        @showtime global sol = sim(u0, p, fdtd_configs)

        # save movies of simulation
        recordsim(sol, p, fdtd_configs, "$name-post_$nres.gif", title="after some training"; frameat, framerate)
    end
end



# plotmonitors(sol0, monitor_instances,)
