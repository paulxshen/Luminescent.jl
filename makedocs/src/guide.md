# Guide
## Implementation
 Length and time are in units of wavelength and period. This normalization allows usage of relative  permitivity and permeability  in equations . Fields including electric, magnetic and current density are simply bundled as a vector of vectors of arrays . Boundary conditions pad the field arrays . PML paddings are multilayered, while All other boundaries add single layers. Paddings are stateful and permanent, increasing the size of field and geometry arrays.  Finite differencing happens every update maxwell_update and are coordinated to implictly implement a staggered Yee's grid .

## Sources
If a source has fewer nonzero dimensions than the simulation domain, its signal will get normalized along its singleton dimensions. For example, all planar sources in 3d or line sources in 2d will get scaled up by a factor of `1/dx`. This way, discretisation would not affect radiated power_flux.
```@docs
PlaneWave
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
power_flux
power_flux_density
```

 ## Physics 
```@docs
maxwell_update!
maxwell_update
```
## GPU support 
Simply use `Flux.gpu` to move simulation variables to GPU. This turns `Arrays` into CUDA arrays which get processed on GPU for both forward and backpropagation passes
## Automatic differentiation adjoints
Compatible with `Zygote.jl` and `Flux.jl`. Please use `maxwell_update` instead of `maxwell_update!` when doing AD. See inverse design examples 
## Generative inverse design
Please contact us for latest scripts. We wrote [Jello.jl](https://github.com/paulxshen/Jello.jl), an innovative Fourier domain neural model to generate length scale controlled efficiently  paramaterized geometry .
## Comparison with other FDTD software
Our focus is on inverse design and topology optimization using adjoints from automatic differentiation. We love a streamlined API backed by a clear, concise and extensible codebase. Attention is paid to speed and malloc but that's not our main concern.
Meep is the most popular open source  FDTD package owing to its maturity and comprehensive features. It supports AD in certain cases using custom adjoints which we avoid in our package in favor of more flexibility .
There are numerous commercial software eg Ansys Lumerical, Comsol and various EDA or SI software. TMK none of them supports AD for inverse design 
