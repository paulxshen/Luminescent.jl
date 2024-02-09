# FDTDEngine.jl
Only 3d works in latest patch. 2d/1d will be fixed in future. Prerelease. Expect breaking changes
## Overview
Differentiable FDTD package for inverse design & topology optimization in photonics, acoustics and RF. Uses automatic differentiation by `Zygote.jl` for adjoint optimization. Integrates with `Jello.jl` to generate length scale controlled paramaterized geometry . Staggered Yee grid update with fully featured boundary conditions & sources. Customizable physics to potentially incorporate dynamics like heat transfer, charge transport.
## Gallery
### Periodic scattering
![](assets/3d_scattering_nres_16.mp4)
### Quarter wavelength antenna
![](assets/3d_quarter_wavelength_antenna_nres_16.mp4)
### Inverse design of compact silicon photonic splitter (coming soon)
In progress, split ratio isn't correct
![](assets/pre_training_nres_16.mp4)
![](assets/post_training_nres_16.mp4)

## Quickstart
We do a quick 3d simulation of plane wave scattering on periodic array of dielectric spheres (first gallery movie)
```julia
using UnPack, LinearAlgebra, GLMakie
using FDTDEngine,FDTDToolkit

F = Float32
name = "3d_scattering"
T = 8.0f0 # simulation duration in [periods]
nres = 16
dx = 1.0f0 / nres # pixel resolution in [wavelengths]
Courant = 0.8 / √3 # Courant number

"geometry"
l = 2 # domain physical size length
sz = nres .* (l, l, l) # domain voxel dimensions
ϵ1 = 1 #
ϵ2 = 2.25f0 # 
b = F.([norm(v .- sz ./ 2) < 0.5 / dx for v = Base.product(Base.oneto.(sz)...)]) # sphere
ϵ = ϵ2 * b + ϵ1 * (1 .- b)

"setup"
boundaries = [Periodic(2), Periodic(3)]# unspecified boundaries default to PML
# n = [1, 0, 0]
monitors = []
sources = [
    PlaneWave(t -> cos(F(2π) * t), -1; Jz=1) # Jz excited plane wave from -x plane (eg -1)
    # PlaneWave(t -> t < 1 ? cos(F(2π) * t) : 0.0f0, -1; Jz=1)
]
configs = setup(boundaries, sources, monitors, dx, sz; F, Courant, T)
@unpack μ, σ, σm, dt, geometry_padding, geometry_splits, field_padding, source_effects, monitor_instances, u0, power = configs

ϵ, μ, σ, σm = apply(geometry_padding; ϵ, μ, σ, σm)
p = apply(geometry_splits; ϵ, μ, σ, σm)


# run simulation
t = 0:dt:T
sol = similar([u0], length(t))
@showtime sol = accumulate!((u, t) -> step!(u, p, t, dx, dt,field_padding, source_effects), sol, t, init=u0)

# make movie
Ez = map(sol) do u
    u[3]
end
ϵz = p[1][3]
dir = @__DIR__
° = π / 180
recordsim(Ez, ϵz, configs, "$dir/$(name)_nres_$nres.mp4", title="$name"; elevation=30°, playback=1, bipolar=true)

```
## Installation
Install via 
```
Pkg.add(url="https://github.com/paulxshen/FDTDEngine.jl")
Pkg.add(url="https://github.com/paulxshen/FDTDToolkit.jl")
```
`FDTDToolkit.jl` contains visualization utilities

## Implementation
Supports 1d (Ez, Hy), 2d TMz (Ez, Hx, Hy), 2d TEz (Hz, Ex, Ey) and 3d. Length and time are in units of wavelength and period. This normalization allows usage of relative  permitivity and permeability  in equations . Fields including electric, magnetic and current density are simply bundled as a vector of arrays . Boundary conditions pad the field arrays . PML paddings are multilayered, stateful and permanent, increasing size of field and geometry arrays. All other boundaries only add transient single layers which are subsequently consumed by finite differencing  every update step. Paddings are coordinated to implictly implement staggered Yee's grid for finite differencing.

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
```

 ## Physics 
```@docs
step3!
```
<!-- step1 -->
<!-- stepTMz -->
<!-- stepTEz -->
## Automatic differentiation adjoints
Pass `bufferfrom=Zygote.bufferfrom` to `step!`
## Tutorials
see `examples/`
## Generative inverse design
Please contact us for latest scripts.
## Community
Discussion & updates at [Julia Discourse](https://discourse.julialang.org/t/pre-ann-differentiable-fdtd-for-inverse-design-in-photonics-acoustics-and-rf/105405/12)
## Contributors
Paul Shen <pxshen@alumni.stanford.edu>
