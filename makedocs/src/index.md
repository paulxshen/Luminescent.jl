# Home
 Currently Prerelease. Expect breaking changes. Report bugs on (Github )[ https://github.com/paulxshen/FDTDEngine.jl] - we usually respond within a day
## Overview
Generative design meets Maxwell's Equations. Differentiable FDTD package for inverse design & topology optimization in semiconductor photonics, acoustics and RF. GPU and automatic differentiation (AD) compatible. Uses AD by `Zygote.jl` for adjoint optimization. Integrates with `Jello.jl` to generate length scale controlled paramaterized geometry . Staggered Yee grid update with fully featured boundary conditions & sources. Customizable physics to potentially incorporate dynamics like heat transfer, charge transport.
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
Pkg.add(url="https://github.com/paulxshen/FDTDEngine.jl")
Pkg.add(url="https://github.com/paulxshen/FDTDToolkit.jl")
```
`FDTDToolkit.jl` contains visualization utilities
## Quickstart
We do a quick 3d simulation of plane wave scattering on periodic array of dielectric spheres (see gallery movie)
```julia
"""
simulation of plane wave scattering on periodic array of dielectric spheres
"""

using UnPack, LinearAlgebra, GLMakie
# using FDTDEngine,FDTDToolkit
dir = pwd()
include("$(dir)/src/main.jl")
include("$dir/../FDTDToolkit.jl/src/main.jl")


name = "3d_scattering"
T = 8.0f0 # simulation duration in [periods]
nx = 16
dx = 1.0f0 / nx # pixel resolution in [wavelengths]

"geometry"
l = 2 # domain physical size length
sz = nx .* (l, l, l) # domain voxel dimensions
ϵ1 = ϵmin = 1 #
ϵ2 = 2.25f0 # 
b = F.([norm(v .- sz ./ 2) < 0.5 / dx for v = Base.product(Base.oneto.(sz)...)]) # sphere
ϵ = ϵ2 * b + ϵ1 * (1 .- b)

"setup"
boundaries = [Periodic(2), Periodic(3)]# unspecified boundaries default to PML
sources = [
    PlaneWave(t -> cos(F(2π) * t), -1; Jz=1) # Jz excited plane wave from -x plane (eg -1)
]
n = [1, 0, 0] # normal 
δ = 0.2f0 # margin
# A = (l - δ)^2
lm = 1 # monitor side length
monitors = [
    Monitor([δ, l / 2, l / 2], [0, lm, lm], n,), # (center, dimensions, normal)
    Monitor([l - δ, l / 2, l / 2], [0, lm, lm], n,),
]
configs = setup(boundaries, sources, monitors, dx, sz; ϵmin, T)
@unpack μ, σ, σm, dt, geometry_padding, geometry_splits, field_padding, source_instances, monitor_instances, u0, = configs

ϵ, μ, σ, σm = apply(geometry_padding; ϵ, μ, σ, σm)
p = apply(geometry_splits; ϵ, μ, σ, σm)


# run simulation
t = 0:dt:T
u = [[similar.(a) for a = u0] for t = t]
u[1] = u0
@showtime reduce(
    (u, (u1, t)) -> step3!(u1, u, p, t, dx, dt, field_padding, source_instances),
    zip(u[2:end], t[1:end-1]),
    init=u0)
y = hcat([power.((m,), u) for m = monitor_instances]...)

# make movie
Ez = map(u) do u
    u[1][3]
end
ϵz = p[1][3]
dir = @__DIR__

recordsim("$dir/$(name)_nres_$nx.mp4", Ez, y;
    dt,
    monitor_instances,
    source_instances,
    geometry=ϵz,
    elevation=30°,
    playback=1,
    axis1=(; title="$name\nEz"),
    axis2=(; title="monitor powers"),
)

```
# Guide
## Implementation
<!-- Supports 1d (Ez, Hy), 2d TMz (Ez, Hx, Hy), 2d TEz (Hz, Ex, Ey) and 3d. -->
 Length and time are in units of wavelength and period. This normalization allows usage of relative  permitivity and permeability  in equations . Fields including electric, magnetic and current density are simply bundled as a vector of vectors of arrays . Boundary conditions pad the field arrays . PML paddings are multilayered, while All other boundaries add single layers. Paddings are stateful and permanent, increasing the size of field and geometry arrays.  Finite differencing happens every update step3 and are coordinated to implictly implement a staggered Yee's grid .

## Sources
If a source has fewer nonzero dimensions than the simulation domain, its signal will get normalized along its singleton dimensions. For example, all planar sources in 3d or line sources in 2d will get scaled up by a factor of `1/dx`. This way, discretisation would not affect radiated power.
```@docs
PlaneWave
Source
```
<!-- GaussianBeam -->

## Boundaries
Unspecified boundaries default to PML 
```@docs
Periodic
PML
PEC
PMC
```
## Monitors  
 ```@docs
Monitor
power
power_density
```

 ## Physics 
```@docs
step3!
step3
```
<!-- step1 -->
<!-- stepTMz -->
<!-- stepTEz -->
## GPU support 
Simply use `Flux.gpu` to move simulation variables to GPU. This turns `Arrays` into CUDA arrays which get processed on GPU for both forward and backpropagation passes
## Automatic differentiation adjoints
Compatible with `Zygote.jl` and `Flux.jl`. Please use `step3` instead of `step3!` when doing AD. See inverse design examples 
## Generative inverse design
Please contact us for latest scripts. We wrote (Jello.jl)[https://github.com/paulxshen/Jello.jl], an innovative Fourier domain neural model to generate length scale controlled efficiently  paramaterized geometry .
## Comparison with other FDTD software
Our focus is on inverse design and topology optimization using adjoints from automatic differentiation. We love a streamlined API backed by a clear, concise and extensible codebase. Attention is paid to speed and malloc but that's not our main concern.
Meep is the most popular open source  FDTD package owing to its maturity and comprehensive features. It supports AD in certain cases using custom adjoints which we avoid in our package in favor of more flexibility .
There are numerous commercial software eg Ansys Lumerical, Comsol and various EDA or SI software. TMK none of them supports AD for inverse design 
# Tutorials
see (`examples/`)[https://github.com/paulxshen/FDTDEngine.jl/blob/main/examples/]
# People
## Community
Discussion & updates at [Julia Discourse](https://discourse.julialang.org/t/pre-ann-differentiable-fdtd-for-inverse-design-in-photonics-acoustics-and-rf/105405/12)
## Contributors
Paul Shen <pxshen@alumni.stanford.edu>
Consulting and technical support available 
2024 (c) Paul Shen