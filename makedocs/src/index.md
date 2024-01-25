# FDTDEngine.jl
Only 3d works in latest patch. 2d/1d will be fixed in future. Prerelease. Expect breaking changes
## Overview
Differentiable FDTD package for inverse design & topology optimization in photonics, acoustics and RF. Uses automatic differentiation by `Zygote.jl` for adjoint optimization. Integrates with `Jello.jl` to generate length scale controlled paramaterized geometry . Staggered Yee grid update with fully featured boundary conditions & sources. Customizable physics to potentially incorporate dynamics like heat transfer, charge transport.

## Quickstart
We do a quick 3d simulation of plane wave on periodic array of dielectric spheres
![](assets/3d_scattering_nres_16.mp4)
```julia
"""
simulation of plane wave scattering on periodic array of dielectric spheres
"""

using UnPack, LinearAlgebra, GLMakie
using FDTDEngine
include("$(pwd())/scripts/plot_recipes.jl")


F = Float32
name = "3d_scattering"
T = 4.0f0 # simulation duration in [periods]
nres = 16
dx = 1.0f0 / nres # pixel resolution in [wavelengths]
Courant = 0.25f0 # Courant number

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
@unpack μ, σ, σm, dt, geometry_padding, geometry_splits, field_padding, source_effects, monitor_instances, fields, step, power = configs

ϵ, μ, σ, σm = apply(geometry_padding; ϵ, μ, σ, σm)
p = apply(geometry_splits; ϵ, μ, σ, σm)
u0 = collect(values(fields))

# run simulation
@showtime sol = accumulate((u, t) -> step(u, p, t, configs), 0:dt:T, init=u0)

# make movie
Ez = map(sol) do u
    u[3]
end
ϵz = p[1][3]
dir = @__DIR__
recordsim(Ez, ϵz, configs, "$dir/$(name)_nres_$nres.mp4", title="$name"; playback=1, bipolar=true)

```
<!-- ![m](assets/periodic_array_nres_32.mp4) -->
## Installation
Install via `Pkg.add(url="https://github.com/paulxshen/FDTDEngine.jl")`. You can additionally access plotting and movie making scripts via `include("_your_path/scripts/plot_recipes.jl")` 
## Implementation
Supports 1d (Ez, Hy), 2d TMz (Ez, Hx, Hy), 2d TEz (Hz, Ex, Ey) and 3d. Length and time are in units of wavelength and period. This normalization allows usage of relative  permitivity and permeability  in equations . Fields including electric, magnetic and current density are simply bundled as a vector of arrays . Boundary conditions pad the field arrays . PML paddings are multilayered, stateful and permanent, increasing size of field and geometry arrays. All other boundaries only add transient single layers which are subsequently consumed by finite differencing  every update step. Paddings are coordinated to implictly implement staggered Yee's grid for finite differencing.

## Sources
If a source has fewer nonzero dimensions than the simulation domain, its signal will get normalized along its singleton dimensions. For example, all planar sources in 3d or line sources in 2d will get scaled up by a factor of `1/dx`. This way, discretisation would not affect radiated power.
```@docs
PlaneWave
GaussianBeam
Source
```

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
```

 ## Physics 
```@docs
step1
stepTMz
stepTEz
step3
```

## Tutorials
Hosted on Google Colab
simulation
- [2d_periodic_array](https://colab.research.google.com/drive/1SiP7MvSE4P05uNgostV9WWx3pFQtFYDW?usp=sharing)
inverse_design
- [2d_waveguide_bend](https://colab.research.google.com/drive/1-g6ShK54MbSsAAeE2c5cj5g7OY7CXfpy?usp=sharing)
## TODO
- integration with external mode solver (suggestions welcome ) for modal source profiles 
- far field transforms
- your wonderful suggestion ;)

## Community
Discussion & updates at [Julia Discourse](https://discourse.julialang.org/t/pre-ann-differentiable-fdtd-for-inverse-design-in-photonics-acoustics-and-rf/105405/12)
## Contributors
Paul Shen <pxshen@alumni.stanford.edu>
