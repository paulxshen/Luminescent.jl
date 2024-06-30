# Guide
## Implementation
 Supports 1d (Ez, Hy), 2d TE (Hz, Ex, Ey), 2d TM (Ez, Hx, Hy), and 3d. Length and time are in units of wavelength and period. This normalization allows usage of relative  permitivity and permeability  in equations . Fields including electric, magnetic and current density are simply bundled as a vector of vectors of arrays . Boundary conditions pad the field arrays . PML paddings are multilayered, while All other boundaries add single layers. Paddings are stateful and permanent, increasing the size of field and geometry arrays.  Finite differencing happens every update update and are coordinated to implictly implement a staggered Yee's grid .

## Sources
We support plane waves and custom  (eg modal) sources . A source has a time signal and a spatial profile . Both can be complex valued in which case the real component is taken from product of the two. Set excited current fields and their spatial profiles as keywords . A profile can be constant , a spatial function , or an array . In the case of array, it'll get resized automatically to fit the source 's spatial extent (spanning from lower bounds `lb` to upper bounds `ub`). 

If a source has fewer nonzero dimensions (from `ub - lb`) than the simulation domain, its signal will get normalized along its singleton dimensions. For example, all planar sources in 3d will get scaled up by a factor of `1/dx`. This way, discretisation would not affect radiated power.

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
field
power
flux
```

 ## Physics 
```@docs
update!
update
```
## GPU support 
Simply use `Flux.gpu` to move simulation variables to GPU. This turns `Arrays` into CUDA arrays which get processed on GPU for both forward and backpropagation passes
## Automatic differentiation adjoints
Compatible with `Zygote.jl` and `Flux.jl`. Please use `update` instead of `update!` when doing AD. See inverse design examples 
## Generative inverse design
Please contact us for latest scripts. We wrote [Jello.jl](https://github.com/paulxshen/Jello.jl), an innovative neural model to generate length scale controlled, efficiently  paramaterized geometry .
## Comparison with other FDTD software
Our focus is on inverse design and topology optimization using adjoints from automatic differentiation. We love a streamlined API backed by a clear, concise and extensible codebase. Attention is paid to speed and malloc but that's not our main concern.
Meep is the most popular open source  FDTD package owing to its maturity and comprehensive features. It supports AD in certain cases using custom adjoints which we avoid in our package in favor of more flexibility .
There are numerous commercial software eg Ansys Lumerical, Comsol and various EDA or SI software. TMK none of them supports AD for inverse design 
