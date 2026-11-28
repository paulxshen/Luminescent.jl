# SOI directional coupler: Tidy3D comparison

Replicates Tidy3D's [directional coupler](https://www.flexcompute.com/tidy3d/examples/notebooks/DirectionalCoupler/) example. Silicon arms (n = 3.48, 0.22 µm thick, 30° sidewalls) on oxide (n = 1.45), air above. Layout from Tidy3D's `coupler.gds`, dilated by 20 nm per side. Light enters port 1 (lower left). We find the wavelength where power splits evenly between ports 3 and 4. Tidy3D reports 1.53 µm with total transmission 0.990 at equivalent nres=15.

```python
import os
import numpy as np
import matplotlib.pyplot as plt

import gdsfactory as gf
from gdsfactory.technology import LogicalLayer, LayerLevel, LayerStack
import luminescent as lumi
lumi.app()
here = os.path.dirname(os.path.abspath(__file__))

# Tidy3D's coupler.gds: 0.45 um arms, 0.1 um gap, 4 um coupling length
gds_file = os.path.join(here, "coupler.gds")
dilation = 0.02  # per side, at waveguide bottom
sidewall_angle = 30  # deg from vertical
h_wg = 0.22

# domain, port positions [um]
lx, ly = 17, 7
x_port = 8
y_port = 2
source_offset = 1

# 21 frequencies over 1.4 - 1.7 um
wavelengths = np.sort(1 / np.linspace(1 / 1.7, 1 / 1.4, 21)).tolist()
wavelength = 2 / (1 / 1.4 + 1 / 1.7)

def setup(path,  sidewall_angle=sidewall_angle):
    nres=8
    shutoff = 1e-3
    tmax = 1000  # max time [periods]
    PML_decay = [0.001, 0.01, 0.1]
    shot_every = 10
    boundaries = ["PML", "PML", "PML"]
    zmin, zmax = -1.11, 1.11
    views = [lumi.view("Hz", z=h_wg / 2)]

    background = "air"

    WG = 1, 0
    BOX = 2, 0
    layer_stack = LayerStack(
        layers={
            "core": LayerLevel(
                layer=LogicalLayer(layer=WG),
                thickness=h_wg,
                zmin=0,
                material="Si",
                mesh_order=0,
                sidewall_angle=sidewall_angle,
            ),
            "box": LayerLevel(
                layer=LogicalLayer(layer=BOX),
                thickness=-zmin,
                zmin=zmin,
                material="SiO2",
                mesh_order=1,
            ),
        }
    )

    # import arms, dilate, crop to domain (extended left to fit source)
    gds = gf.import_gds(gds_file)
    dbu = gds.kcl.dbu
    arms = gf.kdb.Region(gds.begin_shapes_rec(gds.kcl.layer(1, 0)))
    arms = arms.sized(round(dilation / dbu))
    xmin = -max(lx / 2, x_port + source_offset)
    window = gf.kdb.Region(gf.kdb.Box(round(xmin / dbu), round(-ly / 2 / dbu), round(lx / 2 / dbu), round(ly / 2 / dbu)))
    c = gf.Component()
    c.add_polygon(arms & window, layer=WG)
    c.add_polygon(window, layer=BOX)

    # o1: input (lower left), o2: upper left, o3: upper right, o4: lower right
    w_wg = 0.45 + 2 * dilation
    for name, x, y, a in [
        ("o1", -x_port, -y_port, 180),
        ("o2", -x_port, y_port, 180),
        ("o3", x_port, y_port, 0),
        ("o4", x_port, -y_port, 0),
    ]:
        c.add_port(name, center=(x, y), width=w_wg, orientation=a, layer=WG)
    c.plot(return_fig=True).savefig(os.path.join(here, "layout.png"))

    bbox = c.bbox_np().tolist()
    bbox[0].append(zmin)
    bbox[1].append(zmax)

    frames = [
        lumi.plane_frame(
            "wg",
            start=[-1, -0.5],
            stop=[1, zmax],
            wls_mode_solver=[1.4, wavelength, 1.7],
        )
    ]
    sources = [
        lumi.gf_source(
            "o1", c, layer_stack, source_offset,
            frame="wg",
            wl=wavelength,
            bandwidth=wavelength / 8,
        )
    ]
    monitors = [lumi.gf_monitor(p, c, layer_stack, frame="wg") for p in ["o1", "o2", "o3", "o4"]]

    c.flatten()
    lumi.mesh(path, c, layer_stack, lumi.MATERIALS, zmin, zmax)
    lumi.setup(
        keys=["3,1", "4,1"],  # S-params the GUI shows first: the two outputs
        path=path,
        wl=wavelength,
        wls=wavelengths,
        materials=lumi.MATERIALS,
        background=background,
        boundaries=boundaries,
        bbox=bbox,
        frames=frames,
        sources=sources,
        monitors=monitors,
        nres=nres,
        tmax=tmax,
        shutoff=shutoff,
        PML_decay=PML_decay,
        views=views,
        shot_every=shot_every,
    )
    return c
```

```python
path = os.path.join(here, "runs", "SOI_coupler")
setup(path)
```

![layout](img/SOI_coupler_layout.png)

```python
lumi.solve(path)
```

```python
sol = lumi.load(path)
m = sol["monitors"]
a1 = np.array(m["o1"]["modal"][0]["-"])
T3 = abs(np.array(m["o3"]["modal"][0]["+"]) / a1) ** 2
T4 = abs(np.array(m["o4"]["modal"][0]["+"]) / a1) ** 2
wls = np.array(wavelengths)
# 50/50 crossover: linear interpolation where T3 - T4 changes sign
d = T3 - T4
i = np.where(np.diff(np.sign(d)))[0][0]
t = d[i] / (d[i] - d[i + 1])
wl_cross = wls[i] + t * (wls[i + 1] - wls[i])
T_cross = (T3[i] + T4[i]) * (1 - t) + (T3[i + 1] + T4[i + 1]) * t
print(f"crossover at {wl_cross:.3f} um, total transmission {T_cross:.3f}")

plt.plot(wls, T3, label="|S31|²")
plt.plot(wls, T4, label="|S41|²")
plt.axvline(wl_cross, color="gray", ls="--", label=f"crossover {wl_cross:.3f} µm")
plt.xlabel("wavelength (µm)")
plt.ylabel("transmission")
plt.legend()
plt.show()
```

```text
crossover at 1.526 um, total transmission 0.994
```

![out1](img/SOI_coupler_out1.png)

[![GUI](img/SOI_coupler.png)](img/SOI_coupler.png)
