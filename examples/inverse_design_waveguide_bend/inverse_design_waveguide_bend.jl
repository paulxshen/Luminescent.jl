#=
We do inverse design of a compact photonic waveguide bend to demonstrate workflow of FDTD adjoint optimization. First, we seed the design using 2d TE adjoint simulations which serve as fast approximations. Optionlly, we finetune the resulting design in full blown 3d adjoint simulations.

=#

using UnPack, LinearAlgebra, Random, StatsBase, Dates
using Zygote, Flux, CUDA, GLMakie, Jello
using Flux: mae, Adam
using Zygote: withgradient, Buffer
using BSON: @save, @load
using AbbreviatedStackTraces
using Jello, Luminescent, LuminescentVisualization
Random.seed!(1)

# if running directly without module # hide
# include("$(pwd())/src/main.jl") # hide
# include("$(pwd())/../LuminescentVisualization.jl/src/main.jl") # hide

#=
We skip 3d finetuning as it's 20x more compute and memory intensive than 2d adjoints. If wishing to do 3d finetuning, set `iterations3d`. In any case, 3d forward simulations (without adjoint) only take a few seconds.
=#
name = "inverse_design_waveguide_bend"
iterations2d = 10
iterations3d = 0
record2d = true
record3d = false
F = Float32
ongpu = false
model_name = nothing # if load saved model

#=
We load design layout which includes a 2d static_mask of static waveguide geometry as well as variables with locations of ports, signals, design regions and material properties.
=#

@load "$(@__DIR__)/layout.bson" static_mask signals ports designs λ dx ϵsub ϵclad ϵcore hsub hwg hclad
dx, = [dx,] / λ

#=
We initialize a Jello.jl Blob object which will generate geometry of design region. Its parameters will get optimized during adjoint optimization. We initialize it with a straight slab connecting input to output port.
=#

szd = Tuple(round.(Int, designs[1].L / λ / dx) .+ 1) # design region size
if isnothing(model_name)
    nbasis = 5 # complexity of design region
    contrast = 10 # edge sharpness 
    rmin = nothing
    init = [-1 -1 1 -1 -1; -1 1 1 -1 -1; 1 1 -1 -1 -1; -1 -1 -1 -1 -1; -1 -1 -1 -1 -1]
    # init = nothing # random 
    # init = 1 # uniform slab
    model = Blob(szd...; init, nbasis, contrast, rmin,)
else
    @load "$(@__DIR__)/$model_name" model
end
model0 = deepcopy(model)
heatmap(model())

#=
We set key time intervals. The signal must first propagate to port 2 after which all port power fluxes will get monitored
=#

Δ = zeros(2)
# Δ[1] = 1
Δ[1] = 2 + 1.6norm(signals[1].c - ports[2].c) / λ * sqrt(ϵcore) # simulation duration in [periods] for signal to reach output ports
Δ[2] = 2 # duration to record power at output ports
T = cumsum(Δ)

#=
We set boundary conditions, sources , and monitor. The modal source profile is obtained from external mode solver , in our case VectorModesolver.jl . Please refer to guide section of docs website for details . To get an approximate  line source for use in 2d from the cross section profile , we sum and collapse it along its height axis
=#

boundaries = [] # unspecified boundaries default to PML
monitors = [
    # (center, lower bound, upper bound; normal)
    Monitor(p.c / λ, p.lb / λ, p.ub / λ; normal=p.n)
    for p = ports
]

# modal source
@unpack Ex, Ey, Ez, Hx, Hy, Hz = signals[1].modes[1]
Jy, Jx, Mz = map([Ex, Ez, Hy]) do a
    transpose(sum(a, dims=2))
end
Jy, Jx = [Jy, Jx] / maximum(maximum.(abs, [Jy, Jx]))
c = signals[1].c / λ
lb_ = [0, signals[1].lb[1]] / λ
ub_ = [0, signals[1].ub[1]] / λ
sources = [Source(t -> cispi(2t), c, lb_, ub_; Jx, Jy,)]

ϵmin = ϵclad
static_mask = F.(static_mask)
ϵsub, ϵcore, ϵclad = F.((ϵsub, ϵcore, ϵclad))
sz = size(static_mask)

configs = maxwell_setup(boundaries, sources, monitors, dx, sz; F, ϵmin)
@unpack dx, dt, sz, geometry_padding, geometry_staggering, field_padding, source_instances, monitor_instances, u0, = configs

# n = (size(Jy) .- size(monitor_instances[1])) .÷ 2
# power_profile = F.(abs.(Jy[range.(1 .+ n, size(Jy) .- n)...]))
power_profile = F.(real.(Jy .* conj.(Mz)))
power_profile /= norm(power_profile)

if ongpu
    using Flux
    # using CUDA
    # @assert CUDA.functional()
    u0, model, static_mask, μ, σ, σm, field_padding, source_instances =
        gpu.((u0, model, static_mask, μ, σ, σm, field_padding, source_instances))
    merge!(configs, (; u0, field_padding, source_instances))
end

#=
We define a geometry update function that'll be called each adjoint iteration. It calls geometry generator model to generate design region which gets placed onto mask of static features.
    =#

function make_geometry(model, static_mask, configs)#; make3d=false)
    @unpack sz, geometry_padding, geometry_staggering = configs
    μ = ones(F, sz)
    σ = zeros(F, sz)
    σm = zeros(F, sz)
    # μ = 1
    # σ = σm = 0

    mask_ = Zygote.Buffer(static_mask)
    mask_[:, :] = static_mask
    # place!(mask_, σ.(model), round.(Int, designs[1].o / λ / dx) .+ 1)
    place!(mask_, model(), round.(Int, designs[1].o / λ / dx) .+ 1)
    mask = copy(mask_)
    ϵ = mask * ϵcore + (1 .- mask) * ϵclad

    if length(sz) == 3
        ϵ = sandwich(ϵ, round.(Int, [hsub, hwg, hclad] / λ / dx)..., ϵsub, ϵclad)
    end

    p = apply(geometry_padding; ϵ, μ, σ, σm)
    p = apply(geometry_staggering, p)
end

#=
Optimal design will maximize powers into port 1 and out of port 2. Monitor normals were set so both are positive. `metrics` function compute these figures of merit (FOM) quantities by a differentiable FDTD simulation . `loss` is then defined accordingly 
=#

function metrics(model, configs; autodiff=true, history=nothing)
    p = make_geometry(model, static_mask, configs;)
    if !isnothing(history)
        ignore_derivatives() do
            push!(history, p[:ϵ])
        end
    end
    @unpack u0, field_padding, source_instances, monitor_instances = configs
    # run simulation
    _step = if autodiff
        maxwell_update
    else
        maxwell_update!
    end
    u = reduce((u, t) -> _step(u, p, t, dx, dt, field_padding, source_instances;), 0:dt:T[1], init=deepcopy(u0))
    port_fluxes = reduce(T[1]+dt:dt:T[2], init=(u, 0)) do (u, port_fluxes), t
        _step(u, p, t, dx, dt, field_padding, source_instances),
        port_fluxes + dt * flux.((u,), monitor_instances[1:2],)
    end[2] / Δ[2]

    A = support.(monitor_instances)
    port_mode_powers = [mean(vec(a) .* vec(power_profile)) * A for (a, A) = zip(port_fluxes, A)]
    port_powers = mean.(port_fluxes) .* A
    # @info "" port_powers port_mode_powers
    @show port_powers, port_mode_powers
    # println("metrics $port_fluxes")
    abs.(port_mode_powers)
end
# @show const tp = metrics(model, T[1]=1, T[2]=2, autodiff=false)[1] # total power
# error()

function score(v)
    sum(-v)
end

# p0 = make_geometry(model0, static_mask, μ, σ, σm)
history = []
loss = model -> score(metrics(model, configs; history))

#=
We now do adjoint optimization. The first few iterations may show very little change but will pick up momentum
=#

opt = Adam(0.1)
opt_state = Flux.setup(opt, model)
# iterations2d = 66
# iterations2d = 400
for i = 1:iterations2d
    println("$i")
    @time l, (dldm,) = withgradient(loss, model)
    Flux.update!(opt_state, model, dldm)
    println(" $l\n")
end
@save "$(@__DIR__)/2d_model_$(time()).bson" model
# error()

#=
We do a simulation movie using optimized geometry
=#

# @show metrics(model)
function runsave(model, configs; kw...)
    p = make_geometry(model, static_mask, configs)
    @unpack u0, dx, dt, field_padding, source_instances, monitor_instances = configs
    @showtime global u = accumulate((u, t) ->
            maxwell_update!(deepcopy(u), p, t, dx, dt, field_padding, source_instances),
        0:dt:T[2], init=u0)

    # move to cpu for plotting
    if ongpu
        u, p, source_instances = cpu.((u, p, source_instances))
    end
    Hz = field.(u, :Hz)
    ϵEy = field(p, :ϵEy)
    dir = @__DIR__
    d = ndims(Hz[1])
    _name = "$(d)d_$name"
    # error()
    recordsim("$dir/$(_name).mp4", Hz, ;
        dt,
        field=:Hz,
        monitor_instances,
        source_instances,
        geometry=ϵEy,
        rel_lims=0.5,
        playback=1,
        axis1=(; title="$(replace( _name,"_"=>" ")|>titlecase)"),
        axis2=(; title="monitor powers"),
        kw...
    )

end

record = model -> runsave(model, configs)
record2d && record(model)

#=
![](assets/2d_inverse_design_waveguide_bend.mp4)
=#

#=
We now finetune our design in 3d by starting off with optimized model from 2d. We make 3d geometry simply by sandwiching thickened 2d mask between lower substrate and upper clad layers. 
=#

ϵdummy = sandwich(static_mask, round.(Int, [hsub, hwg, hclad] / λ / dx)..., ϵsub, ϵclad)
sz = size(ϵdummy)
model2d = deepcopy(model)


# "monitors"
δ = 0.1 # margin
monitors = [Monitor([p.c / λ..., hsub / λ], [p.lb / λ..., -δ / λ], [p.ub / λ..., hwg / λ + δ / λ]; normal=[p.n..., 0]) for p = ports]

# modal source
@unpack Ex, Ey, Ez, = signals[1].modes[1]
Jy, Jz, Jx = map([Ex, Ey, Ez] / maximum(maximum.(abs, [Ex, Ey, Ez]))) do a
    reshape(a, 1, size(a)...)
end
c = [signals[1].c / λ..., hsub / λ]
lb = [0, signals[1].lb...] / λ
ub = [0, signals[1].ub...] / λ
sources = [Source(t -> cispi(2t), c, lb, ub; Jx, Jy, Jz)]
# sources = [Source(t -> cispi(2t), c, lb, ub; Jx=1)]

configs = maxwell_setup(boundaries, sources, monitors, dx, sz; F, ϵmin, Courant=0.3)
if ongpu
    u0, model, static_mask, μ, σ, σm, field_padding, source_instances =
        gpu.((u0, model, static_mask, μ, σ, σm, field_padding, source_instances))
    merge!(configs, (; u0, field_padding, source_instances))
end


loss = model -> score(metrics(model, configs;))
opt = Adam(0.1)
opt_state = Flux.setup(opt, model)
for i = 1:iterations3d
    @time l, (dldm,) = withgradient(loss, model)
    Flux.update!(opt_state, model, dldm)
    println("$i $l")
end
@save "$(@__DIR__)/3d_model_$(time()).bson" model


record = model -> runsave(model, configs; elevation=70°, azimuth=110°)
record3d && record(model)