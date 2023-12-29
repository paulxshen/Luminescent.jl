# _name_me_please_fdtd.jl
Differentiable FDTD package for inverse design & topology optimization in photonics, acoustics and RF. Uses automatic differentiation by `Zygote.jl` for adjoint optimization. Integrates with `Jello.jl` to generate length scale controlled paramaterized geometry . Staggered Yee grid update with fully featured boundary conditions & sources. Customizable physics to potentially incorporate dynamics like heat transfer, charge transport.

## 

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

## TODO
- integration with external mode solver (suggestions welcome ) for modal source profiles 
- far field transforms
- your wonderful suggestion ;)
 <!-- uses DiffEqFlux.jl neural ODEs and fits within the SciML ecosystem. We have fully featured boundary conditions and support dispersive lossy materials (had to pull some hacks to appease Zygote.jl here). Meant as a differentiable, inverse design focused version of a commercial solver like Ansys Lumerical.
Documentation for _fdtd.jl
I’ll include examples on designing stacks, gratings and couplers as well as meta-materials and meta-surfaces. The electromagnetic FDTD update is simple, drawing on earlier patterns in FDTD.jl and uFDTD.jl based on Schneider’s uFDTD book -->