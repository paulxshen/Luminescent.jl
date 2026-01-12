# Simulation
We recommend first going through a tutorial to get a feel on running simulations. Then come back here for detailed explanations.
## Overview
Everything revolves around the simulation folder `path` you specify in `make`, which uses Python frontend to write `problem.json` and other config files. `solve` then reads them and invokes compiled Julia backend to run simulation, producing `solution.json`. `load` then loads it into a 
## Units
Length and frequency units are arbitrary so long as they are consistent. Time is in periods of the characteristic wavelength. Under the hood everything gets normalized around the characteristic wavelength and its period. 
## Geometry
### Option 1: `gdsfactory`
`gdsfactory` is a popular Python library for programmatic electronic layout. It uses `layer_stack` to convert 2D layers to 3D geometry. `make` accepts `gdsfactory` `layer_stack` and `component` whose ports automatically become simulation ports. .gds and .gerber files can also be imported into `gdsfactory` components.
### Option 2: 3D mesh file import
Create geometry in your favorite CAD and or another simulation tool. Then export each body as .obj (can be saved from .stl or .step) file. Then place them in `geometry` folder inside simulation folder using the convention `{mesh_order}{SEP}{material_name}{SEP}{layer_name}{SEP}{body_name}.obj` eg `1#Si#core#wg.obj`. Lower mesh order overrides higher mesh order. 
## Boundaries
Default is `["PML", "PML", "PML"]` for x, y, z. Each dimension can be further split into + and - eg `[["PEC", "PML"], ...]`. Boundary options include `"PEC"`, `"PMC"`, `"periodic"`. `relative_pml_depths` can shorten PML depths and reduce computation load along dimensions with small or glazing incident radiation eg `[1, 1, 0.3]` in many circuit setups.
## Mesh
Adaptive graded mesh proportions more points to regions of high refractive index and fine features. It also snaps to body boundaries including thin sheets in RF. `nres` arg in `make` controls number of points per characteristic wavelength adjusted locally for refractive index. `mesh_density` relative to vacuum which is by default the refractive index can also be manually set in material constructors. This is useful for mesh overrides or metal materials with undefined refractive index.
## Materials
`materials_library` dictionary maps material names eg `Si` to their property entries which use constructors below. It must contain the special key `background` denoting the background material. The default `lumi.MATERIALS_LIBRARY` contains common materials and can be modified. Material properties are in relative units so ` ϵ₀ = μ₀ = n₀ = Z₀ = 1`. This means you can simply use relative permittivity. Though for other properties like conductivity you must convert to relative units keeping in mind `377 Ω = 1` and `characteristic wavelength in your length units = 1`.
```{eval-rst}
    .. autofunction:: luminescent.setup.Material
    .. autofunction:: luminescent.setup.PECMaterial
    .. autofunction:: luminescent.setup.PlaceholderMaterial
```
## Ports
Plane ports are automatically extracted from gdsfactory component or imported geometry. gdsfactory derived ports are extended by `lateral_port_margin` and `height_port_margin` args in `make` for creating modal sources and monitors. 

Additional ports (not from gdsfactory or marked in imported geometry) can be created as follows and passed to `make` via `ports`
```{eval-rst}
    .. autofunction:: luminescent.setup.PlanePort
    .. autofunction:: luminescent.setup.SpherePort
``` 
## Modes
```{eval-rst}
    .. autofunction:: luminescent.setup.Mode
```
## Monitors
Each port automatically has a monitor. Plane port gets modal monitor which accumulates DFT fields in time domain to obtain field profiles in frequency domain. Modal decomposition then makes forward and backward wave coefficients of each mode wrt frequency. Sphere radiation port works similarly but doesn't necessarily involve modes. 
## Sources
```{eval-rst}
    .. autofunction:: luminescent.setup.Source
```
## Solution
`solve` calls system command `luminescent` to do simulation. It writes `solution.json` which references numpy arrays in the auxiliary folder `.solution.json`. `load` loads them into a dictionary which `query` can further process into derived quantities eg S-parameters 