# Home
 Currently Prerelease. First stable release planned for mid March. Until then, accuracy not validated. Report bugs on [Github](https://github.com/paulxshen/Luminescent.jl) - we usually respond within a day
## Overview
Generative design meets Maxwell's Equations. Differentiable FDTD package for inverse design & topology optimization in semiconductor photonics, acoustics and RF. GPU and automatic differentiation (AD) compatible. Uses AD by `Zygote.jl` for adjoint optimization. Integrates with [`Jello.jl`](https://github.com/paulxshen/Jello.jl) to generate length scale controlled paramaterized geometry . Staggered Yee grid update with fully featured boundary conditions & sources. Customizable physics to potentially incorporate dynamics like heat transfer, charge transport.
## Gallery
### Quarter wavelength antenna radiating above conductive ground plane
![](assets/quarter_wavelength_antenna.mp4)
### Simulation of coupling into dielectric slab waveguide using modal source 
![](assets/slab_waveguide.mp4)
### Simulation of plane wave scattering on Periodic array
![](assets/periodic_scattering.mp4)
### Generative Inverse design of compact silicon photonics splitter 
![](assets/inverse_design_signal_splitter.mp4)

## Installation
Install via 
```
Pkg.add(url="https://github.com/paulxshen/Luminescent.jl")
Pkg.add(url="https://github.com/paulxshen/LuminescentVisualization.jl")
```
`LuminescentVisualization.jl` contains visualization utilities
## Quickstart
We do a quick 3d simulation of plane wave scattering on periodic array of dielectric spheres (see gallery movie)
```julia
"""
simulation of plane wave scattering on periodic array of dielectric spheres
"""

using UnPack, LinearAlgebra, GLMakie
using Luminescent,LuminescentVisualization
dogpu = true
# dogpu = false

name = "periodic_scattering"
T = 10 # simulation duration in [periods]
nx = 20
dx = 1.0 / nx # pixel resolution in [wavelengths]

# geometry
l = 2 # domain physical size length in [wavelengths]
sz = nx .* (l, l, l) # domain voxel dimensions
ϵ1 = ϵmin = 1 #
ϵ2 = 2.25 # 
b = F.([norm(v .- sz ./ 2) < 0.5 / dx for v = Base.product(Base.oneto.(sz)...)]) # sphere
ϵ = ϵ2 * b + ϵ1 * (1 .- b)

# maxwell_setup
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
configs = maxwell_setup(boundaries, sources, monitors, dx, sz; ϵmin, T)
@unpack μ, σ, σm, dt, geometry_padding, geometry_staggering, field_padding, source_instances, monitor_instances, u0, = configs

ϵ, μ, σ, σm = apply(geometry_padding; ϵ, μ, σ, σm)
p = apply(geometry_staggering; ϵ, μ, σ, σm)

# move to gpu
if dogpu
    using CUDA, Flux
    @assert CUDA.functional()
    u0, p, field_padding, source_instances = gpu.((u0, p, field_padding, source_instances))
end

# run simulation
@showtime u = accumulate(0:dt:T, init=u0) do u, t
    maxwell_update!(deepcopy(u), p, t, dx, dt, field_padding, source_instances)
end
y = [power_flux.((m,), u) for m = monitor_instances]

# move back to cpu for plotting
if dogpu
    u, p, field_padding, source_instances = cpu.((u, p, field_padding, source_instances))
end

# make movie, 
Ez = map(u) do u
    u[1][3]
end
ϵEz = p[1][3]
dir = @__DIR__
recordsim("$dir/$(name).mp4", Ez, y;
    dt,
    field=:Ez,
    monitor_instances,
    source_instances,
    geometry=ϵEz,
    elevation=30°,
    playback=1,
    axis1=(; title="$name\nEz"),
    axis2=(; title="monitor powers"),
)
```
