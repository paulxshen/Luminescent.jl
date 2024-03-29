# Home

 Prerelease. First stable release planned for April . Until then, accuracy not validated. Report bugs on [Github](https://github.com/paulxshen/Luminescent.jl) - we usually respond within a day
## Overview
Generative design meets Maxwell's Equations. Differentiable FDTD package for inverse design & topology optimization in semiconductor photonics, acoustics and RF. GPU and automatic differentiation (AD) compatible. Uses AD by `Zygote.jl` for adjoint optimization. Integrates with [`Jello.jl`](https://github.com/paulxshen/Jello.jl) to generate length scale controlled paramaterized geometry . Staggered Yee grid update with fully featured boundary conditions & sources in 1d/2d/3d. Customizable physics to potentially incorporate dynamics like heat transfer, charge transport.
## Gallery

### Generative Inverse design of compact silicon photonics splitter 
![](assets/inverse_design_signal_splitter.mp4)
### Quarter wavelength antenna radiating above conductive ground plane
![](assets/quarter_wavelength_antenna.mp4)
### Simulation of plane wave scattering on Periodic array
![](assets/periodic_scattering.mp4)
Please star us on [Github](https://github.com/paulxshen/Luminescent.jl) if you like our work :)
## Installation
Install via 
```
Pkg.add(url="https://github.com/paulxshen/Luminescent.jl")
Pkg.add(url="https://github.com/paulxshen/LuminescentVisualization.jl")
```
`LuminescentVisualization.jl` contains visualization utilities
## Quickstart
Please refer to the first tutorial 
## People
### Community
Discussion & updates at [Julia Discourse](https://discourse.julialang.org/t/pre-ann-differentiable-fdtd-for-inverse-design-in-photonics-acoustics-and-rf/105405/12)
### Contributors
Paul Shen <pxshen@alumni.stanford.edu>  

Consulting and technical support available  

2024 (c) Paul Shen  