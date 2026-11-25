# Acoplador direccional (TFLN anisótropo): barrido de parámetros



```python
# %%
import os
import numpy as np
import matplotlib.pyplot as plt

import gdsfactory as gf
from gdsfactory.technology import LogicalLayer, LayerLevel, LayerStack
import luminescent as lumi
here = os.path.dirname(os.path.abspath(__file__))

# coupler dims
w_wg = 1
gap = 0.25
# length: to be swept 
dx = 15
dy = 5

wavelength = 1.55
bandwidth = 0.1
nwls = 11
wavelengths = np.linspace(wavelength - bandwidth / 2, wavelength + bandwidth / 2, nwls)

def setup(path, wavelength, wavelengths,  length,gap, dx, dy, w_wg):

    nres = 6
    shutoff = 3e-3
    tmax = 1000  # max time [periods]
    PML_decay = [0.001, 0.01, 0.1]
    shot_every = 10

    h_rib = 0.4
    boundaries = ["PML", "PML", "PML"]
    views = [lumi.view("Hz", x=0, z=h_rib / 2)]

    # layer stack dims
    h_slab = 0.2
    margin = 1.5
    source_offset = 2

    sidewall_angle = 15
    zmargin = 1.
    zmin = -h_slab - zmargin
    zmax = h_rib + zmargin

    # anisotropic TFLN: X-cut along y
    ne = 2.14
    no = 2.21
    materials = {
        "LN": lumi.material(epsilon=[[no**2, 0, 0], [0, ne**2, 0], [0, 0, no**2]]),
        "SiO2": lumi.material(1.44),
    }
    background = "SiO2"

    # layers
    WG = 1, 0
    BBOX = 2, 0
    layer_stack = LayerStack(
        layers={
            "rib": LayerLevel(
                layer=LogicalLayer(layer=WG),
                thickness=h_rib,
                zmin=0,
                material="LN",
                mesh_order=0,
                sidewall_angle=sidewall_angle,
            ),
            "slab": LayerLevel(
                layer=LogicalLayer(layer=BBOX),
                thickness=h_slab,
                zmin=-h_slab,
                material="LN",
                mesh_order=1,
            ),
        }
    )

    # make component in gdsfactory. alternatively can import .gds into gdsfactory
    c = gf.Component()

    dut = c << gf.components.coupler(
        gap=gap,
        length=length,
        dx=dx,
        dy=dy,
        cross_section=gf.cross_section.cross_section(width=w_wg),
    )
    s = c << gf.components.straight(length=source_offset, width=w_wg)
    s.connect("o2", dut.ports["o1"])

    c << gf.components.bbox(component=c, layer=BBOX, top=margin, bottom=margin)

    for p in dut.ports:
        c.add_port(p.name, port=p)
    c.flatten()
    c.plot(return_fig=True).savefig(os.path.join(here, "layout.png"))

    bbox = c.bbox_np().tolist()
    bbox[0].append(zmin)
    bbox[1].append(zmax)

    lx = w_wg + 2 * margin
    frames = [
        lumi.plane_frame(
            "wg",
            start=[-lx / 2, zmin],
            stop=[lx / 2, zmax],
            wls_mode_solver=[wavelength],
        )
    ]
    sources = [
        lumi.gf_source(
            "o1", c, layer_stack, source_offset,
            frame="wg",
            wl=wavelength,
            bandwidth=max(max(wavelengths) - min(wavelengths), wavelength / 8),
        )
    ]
    monitors = [
        lumi.gf_monitor("o1", c, layer_stack, frame="wg"),
        lumi.gf_monitor("o3", c, layer_stack, frame="wg"),
        lumi.gf_monitor("o4", c, layer_stack, frame="wg"),
    ]

    lumi.mesh(path, c, layer_stack, materials, zmin, zmax)
    lumi.setup(
        path=path,
        wl=wavelength,
        wls=wavelengths,
        materials=materials,
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

# %%
S31abs2_length = []
S41abs2_length = []
lengths = [7,8,9]
for length in lengths:
    path = os.path.join(here, "runs", f"TFLN_coupler_{length}")
    setup(path, wavelength, wavelengths, length,gap, dx, dy, w_wg)
    lumi.solve(path)
    sol = lumi.load(path)
    o1 = sol["monitors"]["o1"]
    o3 = sol["monitors"]["o3"]
    o4 = sol["monitors"]["o4"]

    mid=nwls//2
    S31 = o3["modal"][0]["+"] / o1["modal"][0]["-"]
    S41 = o4["modal"][0]["+"] / o1["modal"][0]["-"]

    S31abs2 = np.abs(S31[mid])**2
    S41abs2 = np.abs(S41[mid])**2

    print(f"length: {length} um, S31 power 1.55um: {S31abs2}, S41 power 1.55um: {S41abs2}")
    S31abs2_length.append(S31abs2)
    S41abs2_length.append(S41abs2)
```

![layout](img/TFLN_coupler_layout.png)

