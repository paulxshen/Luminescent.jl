# Simulation Guide
## Workflow
1. Draw the device as a gdsfactory component (or import a .gds) with ports, plus a `LayerStack` assigning thickness and material per layer.
2. Define `materials`, `frames`, `sources` and `monitors`.
3. `lumi.mesh(path, c, layer_stack, materials, zmin, zmax)` writes 3D meshes into simulation folder `path`.
4. `lumi.setup(path, ...)` writes the simulation problem and solves port modes.
5. `lumi.solve(path)` runs FDTD. Runs on Lumi Hub server if `LUMI_API` is set.
6. `lumi.movie(path)` shows the field movie; `sol = lumi.load(path)` loads results.

See [bend](bend.ipynb) for a minimal example.

## Conventions
### Units & frequency
- Length units are arbitrary but must be consistent: eg um for photonics, mm for RF.
- Photonics: `wl` (center) and `wls` (result wavelengths) in length units.
- RF: `f` and `fs` with `c0`, speed of light in your units, eg `c0=300` for mm and GHz.
- Results are arrays over `wls` or `fs`, in the order given.
- Times (`tmax`, `shot_every`, `t0`) are in periods of the center frequency.
- One run computes the whole spectrum; more wavelengths add little runtime.

### Material properties
- `lumi.material(n)` or `lumi.material(epsilon=...)`: relative permittivity, scalar or 3x3 tensor for anisotropic materials.
- Also `mu`, `sigma` (relative units: multiply S/m by `377 * (characteristic wavelength in m)`), `chi2`, `chi3`, `dispersion` (Meep convention).
- `n` also sets mesh density relative to vacuum, default `sqrt(epsilon * mu)`. Raise it to refine the grid in a material, eg thin metals.
- `lumi.PEC_material(n)`: perfect conductor for metals.
- `lumi.placeholder_material(n)`: meshed but no physical properties, eg design regions.
- `materials` is a dict of name to material; `LayerStack` layers and `background` refer to these names.

### Grid & resolution
- `nres`: grid points per wavelength in material. 6 is a good start for photonics; RF usually needs more to resolve thin features.
- Grid adapts to each material's `n`; `subpixel_smoothing` averages material boundaries for accuracy.
- Check convergence by increasing `nres` until results stop changing.
- Use symmetry to halve the domain: PMC or PEC boundary on one side, eg `("PMC", "PML")`.
- Leave `dtype` unset; it's chosen automatically.

### Frames
- Coordinate frame with sampling regions or points, reusable by any number of sources and monitors at different positions via `frame="name"`. Pass all frames to `setup(frames=...)`.
- Local coordinates: origin at the source or monitor `origin`, axes from its basis. `mirrored=True` flips local x, eg for a port on the opposite side of a symmetric device.
- Types and results:
  - `plane_frame`: port cross section with solved modes. Results `modes`, `flux`, `power`. Sources need one.
  - `points_frame`: fields at points, or `sensor="voltage"` / `"current"` integrals for RF.
  - `box_frame`: fields in a sampled volume.
  - `cage_frame`: closed box surface. Flux per face, power, far field.
  - `spherical_frame`: sphere of radius `r`. Flux, power, far field if closed.
- Results are frequency domain at `wls` or `fs` by default; `time_domain=True` records time signals instead.

### Ports & S-parameters
- Source and monitor names match gdsfactory port names, eg `"o1"`. `gf_source` and `gf_monitor` place them at ports.
- General constructors `source` `monitor`
- Port basis rows: local x, local y, outward normal. Sources inject bidirectional.
- Monitor results: `sol["monitors"][name]["modes"][mode number]["+" or "-"]`, where `"+"` is outgoing (along normal), `"-"` incoming.
- S-parameters are ratios of mode amplitudes, eg
```python
S21 = sol["monitors"]["o2"]["modes"][0]["+"] / sol["monitors"]["o1"]["modes"][0]["-"]
S11 = sol["monitors"]["o1"]["modes"][0]["+"] / sol["monitors"]["o1"]["modes"][0]["-"]
```
- Inverse design targets use `"port@mode±"`, eg `lumi.target("o1@0-", "o2@0+", 1)`.

### Inverse design
Optimizes the material layout in design regions to meet S-parameter targets. Not in the free tier. See [GC_invd](GC_invd.ipynb).
1. Draw the design region on its own layer in the component.
2. Give that layer the material `lumi.placeholder_material(n)` in the `LayerStack`, with the lowest `mesh_order` so it overrides other layers, and the vertical extent to be patterned.
3. `lumi.design(name, layer, lmin, fill_material, void_material, ...)`: region patterned with `fill_material` and `void_material` at minimum feature size `lmin`. `repdims=[3]` extrudes a planar pattern along z; `symdims` enforces symmetry.
4. `lumi.target(src, dst, target, wl=...)`: one per objective, eg `"o2@0-"` to `"o1@0+"` with `target=1` for full transmission. Loss is the weighted sum of target errors; use one target per wavelength for broadband.
5. `lumi.optimizer(iters, stoploss, ...)`: stops after `iters` or below `stoploss`. `contrast` below 1 allows topology change and rises quadratically to 1 (binary) as loss approaches `lowloss`, or by iteration 25.
6. Pass `designs`, `targets`, `optimizer` to `setup`, then `solve` runs the optimization.

Outputs in the run folder: `history.json` (loss per iteration) and `designs/<name>/mask.npy` (optimized boolean fill mask).

### Run folder & files
One folder per simulation, eg `runs/bend`:
- `geometry/`: meshes from `mesh`
- `problem.json`, `modes/`, `info.json`, `visualization.json`: from `setup`
- `solution.json`: results from `solve`, read by `load`
- `fields/`, `*.gif`, `*.webp`: field snapshots and movies, if `views` set
