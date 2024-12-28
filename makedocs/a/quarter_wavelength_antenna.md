# Quarter Wavelength Antenna
Complete file at [examples folder](https://github.com/paulxshen/Luminescent.jl/tree/master/examples)


We simulate a quarter wavelength antenna above conductor ground plane and compute its nearfield radiation pattern
```julia

using UnPack, LinearAlgebra, GLMakie, CoordinateTransformations, NearestNeighbors, Random
using GLMakie: volume
using Luminescent, LuminescentVisualization

# if running directly without module # hide
# include("$(pwd())/src/main.jl") # hide
# include("$(pwd())/../LuminescentVisualization.jl/src/main.jl") # hide

Random.seed!
path = "quarter_wavelength_antenna"
F = Float32
dogpu = false
T = 8.0 # simulation duration in [periods]
nx = 30
dx = 1.0 / nx # pixel resolution in [wavelengths]

l = 2 # simulation domain lxlxl box
sz = nx .* (l, l, l)
ϵ = ones(F, sz)
μ = ones(F, sz)
σ = zeros(F, sz)
m = zeros(F, sz)
```
Set Spherical monitor centered on ground. Portions outside domain eg bottom hemisphere are automatically discarded
```julia
boundaries = [PEC(-3)] # ground plane on -z, unspecified boundaries default to PML
monitors = [
    # (center, radius)
    SphereMonitor([l / 2, l / 2, 0], 1),
]
sources = [
    # (signal, center, dimensions)
    Source(t -> cos(2π * t), [l / 2, l / 2, 0.125], [0, 0, 0.25]; Jz=1),
]

prob = setup(boundaries, sources, monitors, dx, sz; F,)
@unpack dt, geometry_padvals, field_lims, field_boundvals, source_instances, monitor_instances, u0, = prob

p = apply(geometry_padvals; ϵ, μ, σ, m)
p = apply(field_lims, p)

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

# move back to cpu for plotting
if dogpu
    u, p, field_boundvals, source_instances = cpu.((u, p, field_boundvals, source_instances))
end
```
Compute nearfield Poynting flux  integrated for 1 period thru our Spherical monitor consisting of points on sphere
```julia
nt = round(Int, 1 / dt)
r = dt * sum(flux.(u[end-nt+1:end], (monitor_instances[1],),))
_, θ, ϕ = eachrow(sphcoords(monitors[1])[:, inbounds(monitor_instances[1])])
```
Interpolate onto regular grid for Plot
```julia
tree = KDTree(hcat(θ, ϕ)')
θ = 0:15°:360°
ϕ = 0:15°:180°
n = length(ϕ)
m = length(θ)
i, = nn(tree, stack(vec(collect.(Base.product(θ, ϕ)))),)
cfs = CartesianFromSpherical()
rvecs = reshape(cfs.(splat(Spherical).(zip(r[i], θ * ones(n)', ones(m) * ϕ'))), m, n)
fig = surface(getindex.(rvecs, 1), getindex.(rvecs, 2), getindex.(rvecs, 3),)
display(fig)
save("$(@__DIR__)/antennapattern.png", fig)
```
![](assets/antennapattern.png)
```julia
```
Ready, set, action! We make movie, 
```julia
Ez = field.(u, :Ez)
ϵEz = field(p, :ϵEz)
dir = @__DIR__
recordsim("$dir/$(name).mp4", Ez, ;
    dt,
    field=:Ez,
    monitor_instances,
    source_instances,
    geometry=ϵEz,
    elevation=30°,
    playback=1,
    axis1=(; title="$name Ez"),
    # axis2=(; title="monitor powers"),
)
```
![](assets/quarter_wavelength_antenna.mp4)
