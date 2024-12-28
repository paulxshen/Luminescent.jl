# Inverse Design Waveguide Bend
Complete file at [examples folder](https://github.com/paulxshen/Luminescent.jl/tree/master/examples)


We do inverse design of a compact photonic waveguide bend to demonstrate workflow of FDTD adjoint optimization. First, we seed the design using 2d TE adjoint simulations which serve as fast approximations. Optionlly, we finetune the resulting design in full blown 3d adjoint simulations.
```julia

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
```
We skip 3d finetuning as it's 20x more compute and memory intensive than 2d adjoints. If wishing to do 3d finetuning, set `iterations3d`. In any case, 3d forward simulations (without adjoint) only take a few seconds.
```julia
path = "inverse_design_waveguide_bend"
iterations2d = 10
iterations3d = 0
record2d = true
record3d = false
F = Float32
ongpu = false
model_path = nothing # if load saved model
```
We load design layout which includes a 2d static_mask of static waveguide geometry as well as variables with locations of ports, sources, design regions and material properties.
```julia

@load "$(@__DIR__)/layout.json" static_mask sources ports designs λ dx ϵbase ϵclad ϵcore hbase hwg hclad
dx, = [dx,] / λ
```
We initialize a Jello.jl Blob object which will generate geometry of design region. Its parameters will get optimized during adjoint optimization. We initialize it with a straight slab connecting input to output port.
```julia

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
```
We set key time intervals. The signal must first propagate to port 2 after which all port power fluxes will get monitored
```julia

deltas = zeros(2)
# deltas[1] = 1
deltas[1] = 2 + 1.6norm(sources[1].c - ports[2].c) / λ * sqrt(ϵcore) # simulation duration in [periods] for signal to reach output ports
deltas[2] = 2 # duration to record power at output ports
T = cumsum(deltas)
```
We set boundary conditions, sources , and monitor. The modal source profile is obtained from external mode solver , in our case VectorModesolver.jl . Please refer to guide section of docs website for details . To get an approximate  line source for use in 2d from the cross section profile , we sum and collapse it along its height axis
```julia

boundaries = [] # unspecified boundaries default to PML
monitors = [
    # (center, lower bound, upper bound; normal)
    Monitor(p.c / λ, p.lb / λ, p.ub / λ; normal=p.n)
    for p = ports
]

# modal source
@unpack Ex, Ey, Ez, Hx, Hy, Hz = sources[1].modes[1]
Jy, Jx, Mz = map([Ex, Ez, Hy]) do a
    transpose(sum(a, dims=2))
end
Jy, Jx = [Jy, Jx] / maximum(maximum.(abs, [Jy, Jx]))
c = sources[1].c / λ
lb_ = [0, sources[1].lb[1]] / λ
ub_ = [0, sources[1].ub[1]] / λ
sources = [Source(t -> cispi(2t), c, lb_, ub_; Jx, Jy,)]

ϵmin = ϵclad
static_mask = F.(static_mask)
ϵbase, ϵcore, ϵclad = F.((ϵbase, ϵcore, ϵclad))
sz = size(static_mask)

prob = setup(boundaries, sources, monitors, dx, sz; F, ϵmin)
@unpack dx, dt, sz, geometry_padvals, field_lims, field_boundvals, source_instances, monitor_instances, u0, = prob

# n = (size(Jy) .- size(monitor_instances[1])) .÷ 2
# power_profile = F.(abs.(Jy[range.(1 .+ n, size(Jy) .- n)...]))
power_profile = F.(real.(Jy .* conj.(Mz)))
power_profile /= norm(power_profile)

if ongpu
    using Flux
    # using CUDA
    # @assert CUDA.functional()
    u0, model, static_mask, μ, σ, m, field_boundvals, source_instances =
        gpu.((u0, model, static_mask, μ, σ, m, field_boundvals, source_instances))
    merge!(prob, (; u0, field_boundvals, source_instances))
end
```
We define a geometry update function that'll be called each adjoint iteration. It calls geometry generator model to generate design region which gets placed onto mask of static features.
    ```julia
function make_geometry(model, static_mask, prob)#; make3d=false)
    @unpack sz, geometry_padvals, field_lims = prob
    μ = ones(F, sz)
    σ = zeros(F, sz)
    m = zeros(F, sz)
    # μ = 1
    # σ = m = 0

    mask_ = Zygote.Buffer(static_mask)
    mask_[:, :] = static_mask
    # place!(mask_, σ.(model), round.(Int, designs[1].o / λ / dx) .+ 1)
    place!(mask_, model(), round.(Int, designs[1].o / λ / dx) .+ 1)
    mask = copy(mask_)
    ϵ = mask * ϵcore + (1 .- mask) * ϵclad

    if length(sz) == 3
        ϵ = sandwich(ϵ, round.(Int, [hbase, hwg, hclad] / λ / dx)..., ϵbase, ϵclad)
    end

    p = apply(geometry_padvals; ϵ, μ, σ, m)
    p = apply(field_lims, p)
end
```
Optimal design will maximize powers into port 1 and out of port 2. Monitor normals were set so both are positive. `metrics` function compute these figures of merit (FOM) quantities by a differentiable FDTD simulation . `loss` is then defined accordingly 
```julia

function metrics(model, prob; autodiff=true, history=nothing)
    p = make_geometry(model, static_mask, prob;)
    if !isnothing(history)
        ignore_derivatives() do
            push!(history, p[:ϵ])
        end
    end
    @unpack u0, field_boundvals, source_instances, monitor_instances = prob
    # run simulation
    _step = if autodiff
        update
    else
        update!
    end
    u = reduce((u, t) -> _step(u, p, t, dx, dt, field_boundvals, source_instances;), 0:dt:T[1], init=deepcopy(u0))
    port_fluxes = reduce(T[1]+dt:dt:T[2], init=(u, 0)) do (u, port_fluxes), t
        _step(u, p, t, dx, dt, field_boundvals, source_instances),
        port_fluxes + dt * flux.((u,), monitor_instances[1:2],)
    end[2] / deltas[2]

    A = area.(monitor_instances)
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

# p0 = make_geometry(model0, static_mask, μ, σ, m)
history = []
loss = model -> score(metrics(model, prob; history))
```
We now do adjoint optimization. The first few iterations may show very little change but will pick up momentum
```julia

opt = RADAM(0.1)
opt_state = Flux.setup(opt, model)
# iterations2d = 66
# iterations2d = 400
for i = 1:iterations2d
    println("$i")
    @time l, (dldm,) = withgradient(loss, model)
    Flux.update!(opt_state, model, dldm)
    println(" $l\n")
end
@save "$(@__DIR__)/2d_model_$(time()).json" model
# error()
```
We do a simulation movie using optimized geometry
```julia

# @show metrics(model)
function runsave(model, prob; kw...)
    p = make_geometry(model, static_mask, prob)
    @unpack u0, dx, dt, field_boundvals, source_instances, monitor_instances = prob
    @showtime global u = accumulate((u, t) ->
            update!(deepcopy(u), p, t, dx, dt, field_boundvals, source_instances),
        0:dt:T[2], init=u0)

    # move to cpu for plotting
    if ongpu
        u, p, source_instances = cpu.((u, p, source_instances))
    end
    Hz = field.(u, :Hz)
    ϵEy = field(p, :ϵEy)
    dir = @__DIR__
    d = ndims(Hz[1])
    _path = "$(d)d_$name"
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

record = model -> runsave(model, prob)
record2d && record(model)
```
![](assets/2d_inverse_design_waveguide_bend.mp4)
```julia
```
We now finetune our design in 3d by starting off with optimized model from 2d. We make 3d geometry simply by sandwiching thickened 2d mask between lower substrate and upper clad layers. 
```julia

ϵdummy = sandwich(static_mask, round.(Int, [hbase, hwg, hclad] / λ / dx)..., ϵbase, ϵclad)
sz = size(ϵdummy)
model2d = deepcopy(model)


# "monitors"
δ = 0.1 # margin
monitors = [Monitor([p.c / λ..., hbase / λ], [p.lb / λ..., -δ / λ], [p.ub / λ..., hwg / λ + δ / λ]; normal=[p.n..., 0]) for p = ports]

# modal source
@unpack Ex, Ey, Ez, = sources[1].modes[1]
Jy, Jz, Jx = map([Ex, Ey, Ez] / maximum(maximum.(abs, [Ex, Ey, Ez]))) do a
    reshape(a, 1, size(a)...)
end
c = [sources[1].c / λ..., hbase / λ]
lb = [0, sources[1].lb...] / λ
ub = [0, sources[1].ub...] / λ
sources = [Source(t -> cispi(2t), c, lb, ub; Jx, Jy, Jz)]
# sources = [Source(t -> cispi(2t), c, lb, ub; Jx=1)]

prob = setup(boundaries, sources, monitors, dx, sz; F, ϵmin, Courant=0.3)
if ongpu
    u0, model, static_mask, μ, σ, m, field_boundvals, source_instances =
        gpu.((u0, model, static_mask, μ, σ, m, field_boundvals, source_instances))
    merge!(prob, (; u0, field_boundvals, source_instances))
end


loss = model -> score(metrics(model, prob;))
opt = RADAM(0.1)
opt_state = Flux.setup(opt, model)
for i = 1:iterations3d
    @time l, (dldm,) = withgradient(loss, model)
    Flux.update!(opt_state, model, dldm)
    println("$i $l")
end
@save "$(@__DIR__)/3d_model_$(time()).json" model


record = model -> runsave(model, prob; elevation=70°, azimuth=110°)
record3d && record(model)
```