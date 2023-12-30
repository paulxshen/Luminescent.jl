# fdtd_prerelease.jl
Prerelease. Expect breaking changes
## Overview
Differentiable FDTD package for inverse design & topology optimization in photonics, acoustics and RF. Uses automatic differentiation by `Zygote.jl` for adjoint optimization. Integrates with `Jello.jl` to generate length scale controlled paramaterized geometry . Staggered Yee grid update with fully featured boundary conditions & sources. Customizable physics to potentially incorporate dynamics like heat transfer, charge transport.

## Quickstart
We do a quick 2d simulation of plane wave on periodic array of dielectric pillars in front of metal of plane
![](assets/periodic_array_nres_32.mp4)
```julia
name = "periodic_array"
using UnPack,fdtd_prerelease
include("examples/utils/plot_recipes.jl")

"simulation params"
T = 20.0f0 # simulation duration in [periods]
nres = 32
Courant = 0.5f0 # Courant number
F = Float32

"recording params"
frameat = 1 / 16 # captures frame every _ periods
framerate = 16 # playback speed

dx = 1.0f0 / nres # pixel resolution in [wavelengths]
L = F[8, 8] # domain dimensions in [wavelengths]

# setup FDTD
polarization = :TMz
boundaries = [Periodic(2), PEC(1)] # unspecified boundaries default to PML
monitors = [Monitor([L / 2])]
sources = [PlaneWave(t -> exp(-((t - 1) / 1)^2) * cos(F(2π) * t), -1; Jz=1)]
fdtd_configs = setup(boundaries, sources, monitors, L, dx, polarization; F, Courant, T)
@unpack esz0, hsz0, dt, geometry_padding, field_padding, source_effects, monitor_instances, fields = fdtd_configs

ϵ1 = 1 #
ϵ2 = 2.25f0 # 
b = F.([norm([x, y] .- esz0 ./ 2) < 1 / dx for x = 1:esz0[1], y = 1:esz0[2]]) # circle

ϵ = ϵ2 * b + ϵ1 * (1 .- b)
μ = ones(F, hsz0)
σ = zeros(F, esz0)
σm = zeros(F, hsz0)
ϵ, μ, σ, σm = apply(geometry_padding; ϵ, μ, σ, σm)

p = [ϵ, μ, σ, σm]
u0 = collect(values(fields))

# run & record simulation
@showtime sol = accumulate((u, t) -> stepTMz(u, p, t, fdtd_configs), 0:dt:T, init=u0)
recordsim(sol, p, fdtd_configs, "$(name)_nres_$nres.mp4", title="$name"; frameat, framerate)
```
<!-- ![m](assets/periodic_array_nres_32.mp4) -->
## Installation
Install via `Pkg.add(url="https://github.com/paulxshen/differentiable-fdtd-beta-prerelease")`. You can additionally access plotting and movie making scripts via `include("examples/utils/plot_recipes.jl")` (may need to modify path)
## Implementation
Supports 1d (Ez, Hy), 2d TMz (Ez, Hx, Hy), 2d TEz (Hz, Ex, Ey) and 3d. Length and time are in units of wavelength and period. This normalization allows usage of relative  permitivity and permeability  in equations . Fields including electric, magnetic and current density are simply bundled as a vector of arrays . Boundary conditions pad the field arrays . PML paddings are multilayered, stateful and permanent, increasing size of field and geometry arrays. All other boundaries only add transient single layers which are subsequently consumed by finite differencing  every update step. Paddings are coordinated to implictly implement staggered Yee's grid for finite differencing.

## Sources
```@docs
PlaneWave
GaussianBeam
UniformSource
CenteredSource
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
