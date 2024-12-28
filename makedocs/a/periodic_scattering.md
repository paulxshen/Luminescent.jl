# Periodic Scattering
Complete file at [examples folder](https://github.com/paulxshen/Luminescent.jl/tree/master/examples)


We simulate  plane wave scattering on periodic array of dielectric spheres
```julia
using UnPack, LinearAlgebra, GLMakie
using Luminescent, LuminescentVisualization

# if running directly without module # hide
# include("$(pwd())/src/main.jl") # hide
# include("$(pwd())/../LuminescentVisualization.jl/src/main.jl") # hide
``` 
Set simulation duration and resolution.  Run on CPU by setting `dogpu = false`. If running on a newer GPU, set `F = Float16`
```julia
path = "periodic_scattering"
T = 10 # simulation duration in [periods]
nx = 20
dx = 1.0 / nx # pixel resolution in [wavelengths]
dogpu = false
F = Float32
```
We make unit cell geometry containing a dielectric sphere. Each property is made an array
```julia
l = 2 # domain physical size length in [wavelengths]
sz = nx .* (l, l, l) # domain voxel dimensions

ϵ1 = ϵmin = 1 #
ϵ2 = 2.25 # 
b = F.([norm(v .- sz ./ 2) < 0.5 / dx for v = Base.product(Base.oneto.(sz)...)]) # sphere
ϵ = ϵ2 * b + ϵ1 * (1 .- b)

# μ = 1
μ = ones(F, sz)
σ = zeros(F, sz)
m = zeros(F, sz)
```
We setup boundary conditions, source and monitor surfaces
```julia
boundaries = [Periodic(2), Periodic(3)]# unspecified boundaries default to PML
sources = [
    PlaneWave(t -> cos(2π * t), -1; Jz=1) # Jz excited plane wave from -x plane (eg -1)
]
normal = [1, 0, 0] #  
δ = 0.2 # margin
lm = 1 # monitor side length
monitors = [
    Monitor([δ, l / 2, l / 2], [0, lm, lm]; normal), # (center, dimensions; normal)
    Monitor([l - δ, l / 2, l / 2], [0, lm, lm]; normal),
]
```
We do `setup` to instantiate at the given discretisation. We adopt `u, p, t` naming conventions from ODE literature: `u ` as state, `p` as params eg geometry
```julia
prob = setup(boundaries, sources, monitors, dx, sz; ϵmin, F)
@unpack dt, geometry_padvals, field_lims, field_boundvals, source_instances, monitor_instances, u0, = prob

p = apply(geometry_padvals; ϵ, μ, σ, m)
p = apply(field_lims; p...)

# move to gpu
if dogpu
    using Flux
    # using CUDA
    # @assert CUDA.functional()
    u0, p, field_boundvals, source_instances = gpu.((u0, p, field_boundvals, source_instances))
end
```
We run simulation as an `accumulate` loop. `update!` applies Maxwells equations as staggered time stepping on E, H. It's mutating so a copy is made in order to save sequence of states
```julia
@showtime u = accumulate(0:dt:T, init=u0) do u, t
    update!(deepcopy(u), p, t, dx, dt, field_boundvals, source_instances)
end
port_powers = [power.(u, (m,),) for m = monitor_instances]

# move back to cpu for plotting
if dogpu
    u, p, field_boundvals, source_instances = cpu.((u, p, field_boundvals, source_instances))
end
```
Ready, set, action! We make movie, 
```julia
Ez = field.(u, :Ez)
ϵEz = field(p, :ϵEz)
dir = @__DIR__
recordsim("$dir/$(name).mp4", Ez, port_powers;
    dt,
    field=:Ez,
    monitor_instances,
    source_instances,
    geometry=ϵEz,
    elevation=30°,
    playback=1,
    axis1=(; title="$name"),
    axis2=(; title="monitor powers"),
)
```
![](assets/periodic_scattering.mp4)
