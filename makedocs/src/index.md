# Home
06/2024 update: we've taken down temporarily the Julia API documentation after significant revamp. Meanwhile , refer to our [Python gdsfactory plugin](luminescent.html) for simulation & inverse design of photonic integrated circuits 
## Overview
Generative design meets Maxwell's Equations. GPU and automatic differentiation (AD) compatible FDTD package in Julia for inverse design & topology optimization in semiconductor photonics, acoustics and RF. Uses AD by `Zygote.jl` for adjoint optimization. Integrates with [`Jello.jl`](https://github.com/paulxshen/Jello.jl) to generate length scale controlled paramaterized geometry . Staggered Yee grid update with fully featured boundary conditions & sources in 1d/2d/3d. 
## Gallery

###  Inverse design of photonic wavelength domain demultiplexer 
![](assets/demux.png)
### Quarter wavelength antenna radiating above conductive ground plane
![](assets/quarter_wavelength_antenna.mp4)
### Simulation of plane wave scattering on Periodic array
![](assets/periodic_scattering.mp4)
Please star us on [Github](https://github.com/paulxshen/Luminescent.jl) if you like our work :)
## People
### Community
Discussion & updates at [Julia Discourse](https://discourse.julialang.org/t/pre-ann-differentiable-fdtd-for-inverse-design-in-photonics-acoustics-and-rf/105405/12)
### Contributors
Paul Shen <pxshen@alumni.stanford.edu>  
Part of Luminescent AI <info@luminescentai.com>

2024 (c) Paul Shen  