# Python API
## Workflow
### mesh
```python
lumi.mesh(path, component, layer_stack, materials, zmin, zmax)
```
Writes .obj meshes of gdsfactory `component` extruded by `layer_stack`, clipped to [zmin, zmax], into simulation folder `path`. Call before `setup`.

### setup
```python
lumi.setup(path, sources, monitors, nres=6, wls=None, wl=None, f=None, fs=None, c0=None, targets=[], designs=[], approx_2D_mode=None, **kwargs)
```
Writes simulation problem to folder `path`. Call `mesh` first, then `solve`.

Args:

- wl and/or wls: center wavelength and wavelengths to compute results at, or
- f and/or fs with c0: center frequency and frequencies, with speed of light c0 in your units eg 300 [mm GHz]. Results are then ordered by frequency
- sources, monitors: from `source`, `monitor`, `gf_source`, `gf_monitor`
- frames: from `plane_frame`, `box_frame`, `spherical_frame`, `points_frame`, `cage_frame`
- materials: dict of name to `material`
- background: material name filling empty space
- bbox: simulation box [[xmin, ymin, zmin], [xmax, ymax, zmax]]
- boundaries: per axis "PML", "PEC", "PMC", "periodic" or pair (min side, max side) eg ("PEC", "PML")
- nres: grid points per wavelength in material. 6 is a good starting point for photonics; increase to check convergence
- tmax: max simulation time [periods]
- stop_decay: stop early once field energy decays to this fraction
- PML_decay: PML field attenuation, scalar or per axis
- relative_courant: time step relative to max stable
- subpixel_smoothing: subpixel averaging of material boundaries
- views: from `view`, for movie
- shot_every: save movie frame every _ periods
- device: "cpu" or "gpu"
- dtype: leave unset; chosen automatically. Pass "f32" or "f16" only to override
- approx_2D_mode: 2D approximation
- designs, targets, optimizer: inverse design, from `design`, `target`, `optimizer`

### solve
```python
lumi.solve(path)
```
Run simulation in folder `path` prepared by `setup`. Writes solution.json.
Runs on a Luminescent server if env LUMI_API is set (see `remote_solve`).
Runs the Julia source instead of the compiled binary if env LUMIDEV is set.

### load
```python
lumi.load(path, show=True)
```
Load solution dict from `path`. Monitor results in sol["monitors"][name], eg ["modes"][mode number]["+" or "-"] as array over wavelengths (or frequencies if setup with fs).

### movie
```python
lumi.movie(path, views=None, fps=None)
```
Show simulation movie in `path`, rendered by `solve`.
`views` (overrides views from setup) or `fps` re-render it from field snapshots, which cloud solves don't keep.

### view
```python
lumi.view(field, name='', prop='n', x=None, y=None, z=None, field_color_intensity=1.0, show_grid=False, show_box=True, labels='auto', mirrors=[])
```
Movie view: slices of `field` overlaid on material property `prop`, for `setup(views=...)`.

Args:

- field: field component eg "Ey", "Hz"
- name: movie file name, defaults to field and slice positions
- prop: material property underneath eg "n", "epsilon"
- x, y, z: slice position along that axis, "half" for the simulation box center, or None for no slice. One panel per slice
- field_color_intensity: scales field color saturation
- show_grid: draw grid lines
- show_box: draw simulation box outline
- labels: "auto" to label sources and monitors, or dict of name to label
- mirrors: reflect across symmetry boundaries for display eg ["z-"]

## Materials
### material
```python
lumi.material(n=None, epsilon=None, mu=1, sigma=0, chi2=0, chi3=0, dispersion=[], placeholder=False)
```
Define material properties in relative units

Args:

- epsilon: relative permittivity
- mu: relative permeability
- sigma: conductivity (multiply by `377 * (characteristic wavelength in m)` to convert from metric (`S/m, 1 /(Ω m)`) to relative units)
- chi2: Pockels second order nonlinear susceptibility
- chi3: Kerr third order nonlinear susceptibility
- dispersion: list of tuples [(gamma, omega, is_lorentzian)...] each defining a dispersion term. `is_lorentzian=True` for Lorentzian, `False` for Drude. We use Meep convention https://meep.readthedocs.io/en/latest/Materials/#material-dispersion and relative units.
- n: mesh density relative to vacuum, defaults to square root of relative permittivity

### PEC_material
```python
lumi.PEC_material(n)
```
Perfect electric conductor material for modeling metals

### placeholder_material
```python
lumi.placeholder_material(n)
```
placeholder material contributes to meshing but not to physical properties, eg used for design region

## Sources & monitors
### source
```python
lumi.source(name, origin, basis, frame, mirrored=False, wl=None, f=None, bandwidth=None, duration=None, modenums=[0])
```
Modulated Gaussian pulse source with gradual spectral roll off  at bandwidth limits. Becomes continuous wave if `bandwidth=0`.

Args:

- name: should match a port name eg "o1"
- origin: position [x, y, z]
- basis: rows local x, local y, outward normal. source injects bidirectionally (along both ±normal)
- frame: name of `plane_frame` holding mode profile
- mirrored: reflect frame and modes across local x=0, eg for port on opposite side of a symmetric frame
- wl or f: wavelength or frequency, specify either (not both). f requires c0 in setup
- bandwidth: in same units as above
- modenums: list of mode numbers to excite

### monitor
```python
lumi.monitor(name, origin, basis, frame=None, mirrored=False, t0=0, x=None, y=None, z=None, direction='+', layer=None)
```
Monitor recording results of `frame` at `origin` in `basis`.

Args:

- name: should match a port name eg "o1". Unique among monitors
- origin: position [x, y, z]
- basis: rows local x, local y, outward normal. "+" results propagate along normal, "-" against
- frame: frame name
- mirrored: reflect frame and modes across local x=0
- t0: start recording time-domain signals at t0 [periods]

### gf_source
```python
lumi.gf_source(port, c, layer_stack, offset=0, name=None, **kwargs)
```
`source` at gdsfactory port `port` of component `c`, `offset` along outward normal. `name` defaults to port name. kwargs go to `source`.

### gf_monitor
```python
lumi.gf_monitor(port, c, layer_stack, offset=0, name=None, **kwargs)
```
`monitor` at gdsfactory port `port` of component `c`. `name` defaults to port name. kwargs go to `monitor`.

### gaussian_mode
```python
lumi.gaussian_mode(r, R, Z, dx, wl, theta=0, half=None)
```
Tilted Gaussian beam mode for `plane_frame(modes=...)`, eg fiber.

Args:

- r: beam radius (1/e field at r*sqrt(2))
- R: frame half width
- Z: wave impedance in vacuum units, eg 1/n
- dx: sampling
- wl: wavelength
- theta: tilt angle [rad] in local xz plane
- half: "y+" for half frame y>=0 (symmetry plane at y=0)

## Frames
### plane_frame
```python
lumi.plane_frame(name, start, stop, fs_mode_solver=[], wls_mode_solver=[], modes='auto', nmodes=1, metallic_boundaries=None, boundaries=('PML', 'PML'), time_domain=None, frequency_domain=None)
```
Plane port frame between local `start`, `stop` (2 coords in 3D). Mode solutions at given frequencies or wavelengths; sources and monitors use (or interpolate) modes at closest wavelengths. Results: modes, flux, power. Records frequency-domain results by default, or time-domain signals if `time_domain=True` (or `frequency_domain=False`).

Args:

- fs_mode_solver or wls_mode_solver: frequencies (requires c0 in setup) or wavelengths to solve modes at, specify one (not both)
- modes: "auto" to solve, else precomputed fields eg [[{'Ex': np_array, 'Hy': ...}, ...], ...]
- nmodes: number of modes to solve for
- metallic_boundaries: PEC outlines for mode solver: layer names or frame edges "x-", "x+", "y-", "y+"

### points_frame
```python
lumi.points_frame(name, points, sensor=None, time_domain=None, frequency_domain=None)
```
Points in monitor basis. Records frequency-domain results by default, or time-domain signals if `time_domain=True` (or `frequency_domain=False`).

Args:

- sensor: "voltage": ∫E·dl along polyline `points`. "current": ∮H·dl around closed loop `points`. None: raw fields.

### box_frame
```python
lumi.box_frame(name, start, stop, dr, time_domain=None, frequency_domain=None)
```
Box volume frame between local `start`, `stop` (3 coords in 3D), sampled every `dr` (scalar or per axis). Results: fields in monitor basis. Records frequency-domain results by default, or time-domain signals if `time_domain=True` (or `frequency_domain=False`).

### cage_frame
```python
lumi.cage_frame(name, start, stop, dr, theta=3.141592653589793, phi=6.283185307179586, dtheta=np.float64(0.08726646259971647), dphi=np.float64(0.08726646259971647), time_domain=None, frequency_domain=None)
```
Cage: closed surface of box between local `start`, `stop` (3D only), faces sampled every `dr` (scalar or per axis). Results: flux per face (outward), power, farfield intensity and directivity on (theta, phi) grid. theta, phi: max angle (range from 0) or [min, max]. Records frequency-domain results by default, or time-domain signals if `time_domain=True` (or `frequency_domain=False`).

### spherical_frame
```python
lumi.spherical_frame(name, r, theta=3.141592653589793, phi=6.283185307179586, dtheta=np.float64(0.08726646259971647), dphi=np.float64(0.08726646259971647), time_domain=None, frequency_domain=None)
```
Sphere (circle in 2D) of radius `r` about monitor origin, in monitor basis. theta, phi: max angle (range from 0) or [min, max]. Results: flux (radial), power, and if closed (full theta, phi), farfield intensity and directivity on (theta, phi) grid. Records frequency-domain results by default, or time-domain signals if `time_domain=True` (or `frequency_domain=False`).

## Inverse design
### design
```python
lumi.design(name, layer, lmin, fill_material, void_material, component, layer_stack, materials, symdims=[], repdims=[3], anchordims=[], dependent=False, topopt=True)
```
Design region to be optimized

Args:

- name: design name eg "d1"
- layer: placeholder layer for design region
- lmin: minimum feature size
- fill_material: eg 'Si'
- void_material: eg 'SiO2'
- symdims: eg ["x","y","diagonal"]
- component, layer_stack, materials: as in setup
- repdims: eg [3] for planar design region extruded along z
- anchordims: signed dims eg [-1, 2] whose design edges (x-, y+) are closed instead of open
- topopt: allow topology change

### target
```python
lumi.target(src, dst, target, target_func='abs2', error_func='abs', weight=1, f=None, wl=None)
```
- src: source port and mode eg "o1@0-" (mode 0 incoming)
- dst: destination port and mode eg "o2@0+" (mode 0 outgoing)
- target: target value
- target_func: applied to S parameter before comparing to target eg "abs2" for power
- weight: weight of this target in overall loss function
- f or wl: frequency or wavelength, specify one (not both)
- error_func: error function `abs, identity, -, relu, relu-'

### optimizer
```python
lumi.optimizer(iters, stoploss, startloss=1, lowloss=None, momentum=0.7, contrast=0.8, init_dA=0.03, ckptat=None, gradckptat=None)
```
Args:

- iters: maximum number of optimization iterations
- stoploss: stop optimization if loss goes below this value
- startloss: initial loss estimate
- lowloss: loss at which contrast reaches 1. defaults to 0.2 * stoploss + 0.8 * startloss
- momentum: optimizer momentum
- contrast: initial binary contrast (0 to 1). slightly <1 enables topology change (new holes or islands). =1 freezes topology eg shape optimization. contrast automatically increases to 1 as loss decreases to `lowloss`
- init_dA: fraction of design area changed by the first step. sets the initial step size
- ckptat: save optimization checkpoint every _ iterations
- gradckptat: gradient checkpointing

### apply_design
```python
lumi.apply_design(c0, sol)
```
Replace design regions of component `c0` with optimized designs from solution `sol`.

### design_from_gds
```python
lumi.design_from_gds(path, i=1)
```
Import `i`th optimized design from `path` as gdsfactory component.

## Utilities
### straight
```python
lumi.straight(length=10, width: float = 1, layer=(1, 0))
```
Rectangular gdsfactory component with ports "o1" (west) and "o2" (east).

### remote_solve
```python
lumi.remote_solve(path, api=None, token=None, poll=3)
```
Run simulation in folder `path` on a Luminescent server. Uploads inputs, streams the log, downloads results into `path`.
The server renders the movie; field snapshots are not downloaded.
`api` and `token` default to env LUMI_API and LUMI_TOKEN (or JUPYTERHUB_API_TOKEN in hosted notebooks).

### load_sparams
```python
lumi.load_sparams(sparams)
```
Convert {"re", "im"} dicts to complex values.

### version
```python
lumi.version()
```
Installed luminescent python package version. Its X.Y must match the solver binary's X.Y.0.

