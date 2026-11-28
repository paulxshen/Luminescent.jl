# Microtira de RF: filtro de stub de cuarto de onda

Filtro de rechazo de banda con stub de cuarto de onda en microtira de RF con muesca nominal a 5GHz

```python
import os
import numpy as np
import matplotlib.pyplot as plt

import gdsfactory as gf
from gdsfactory.technology import LogicalLayer, LayerLevel, LayerStack
import luminescent as lumi
lumi.app()
here = os.path.dirname(os.path.abspath(__file__))

# simulation folder
path = os.path.join(here, "runs", "stub_filter")

# dims [mm]
# length units are arbitrary so long as they're consistent
w_line = 1.52  # strip width
h_metal = 0.035  # metal thickness
h_sub = 0.8  # substrate thickness
v_rel = 87.5 / 139  # relative signal speed, from microstrip calculator

# frequencies [GHz]
# spectrum computed from single FDTD run. more frequencies slightly increases runtime
frequency = 5
bandwidth = 4
nfs = bandwidth * 25+1
frequencies = np.linspace(frequency - bandwidth / 2, frequency + bandwidth / 2, nfs)

c0 = 300  # speed of light [mm GHz]

# quarter-wave stub at center frequency
l_stub = c0 / frequency * v_rel / 4
l_feed = l_stub

def setup(path, frequency, frequencies, w_line, h_metal, h_sub, l_stub, l_feed):
    # simulation parameters
    nres = 8
    dx = h_sub / 6  # desired resolution in metal
    shutoff = 1e-2
    tmax = 1000  # max time [periods]
    PML_decay = [0.001,.1,.1]
    courant = 0.9
    shot_every =.5 # movie [periods]

    # margins
    zmargin =margin = 5 * h_sub
    source_offset = 5 * h_sub

    # ground plane as PEC boundary on z-
    boundaries = ["PML", "PML", ("PEC", "PML")]

    # vertical dims
    zmin = -h_sub
    zmax = h_metal + zmargin

    views = [
        lumi.view("Ez", x=-source_offset),
        lumi.view("Ez", y=0),
        lumi.view("Ez", z=-h_sub / 2),
        lumi.view("Hz", z=h_metal / 2),
    ]

    # materials
    # refine mesh in metal: n sets mesh density relative to vacuum
    n = (c0 / frequency / nres) / dx
    materials = {**lumi.MATERIALS, "PEC": lumi.PEC(n)}
    background = "air"

    # layers
    TRACE = 1, 0
    SUB = 10, 0
    layer_stack = LayerStack(
        layers={
            "trace": LayerLevel(
                layer=LogicalLayer(layer=TRACE),
                thickness=h_metal,
                zmin=0,
                material="PEC",
                mesh_order=10,
            ),
            "sub": LayerLevel(
                layer=LogicalLayer(layer=SUB),
                thickness=h_sub,
                zmin=-h_sub,
                material="FR4",
                mesh_order=30,
            ),
        }
    )

    # make component in gdsfactory. alternatively can import .gds into gdsfactory
    c = gf.Component()
    feed = c << lumi.straight(length=l_feed, width=w_line, layer=TRACE)
    stub = c << lumi.straight(length=l_stub, width=w_line, layer=TRACE)
    stub.rotate(90)
    stub.movex(l_feed / 2)
    stub.movey(-w_line / 2)
    s = c << lumi.straight(length=source_offset, width=w_line, layer=TRACE)
    s.connect("o2", feed.ports["o1"])

    c << gf.components.bbox(component=c, layer=SUB, top=margin, bottom=margin)

    c.add_port("o1", port=feed.ports["o1"])
    c.add_port("o2", port=feed.ports["o2"])
    c.plot(return_fig=True).savefig(os.path.join(here, "layout.png"))

    bbox = c.bbox_np().tolist()
    bbox[0].append(zmin)
    bbox[1].append(zmax)

    # sources, monitors
    # quasi-TEM microstrip mode: ground plane (y-) and strip ("trace" layer) as metal
    w_port = w_line + 2 * margin
    frames = [
        lumi.plane_frame(
            "line_mode",
            start=[-w_port / 2, zmin],  # local xy basis, local y=0 at global zmin of port layer
            stop=[w_port / 2, h_metal + margin],
            metallic_boundaries=["y-", "trace"],
            nmodes=1,
            modes="auto",
            fs_mode_solver=[frequency],
        ),
    ]

    sources = [
        lumi.gf_source(
            "o1", c, layer_stack, source_offset,
            frame="line_mode",
            f=frequency,
            bandwidth=max(max(frequencies) - min(frequencies), frequency / 8),
        )
    ]

    monitors = [
        lumi.gf_monitor("o1", c, layer_stack, frame="line_mode"),
        lumi.gf_monitor("o2", c, layer_stack, frame="line_mode"),
    ]

    # writes .obj mesh files from gdsfactory component
    c.flatten()
    lumi.mesh(path, c, layer_stack, materials, zmin, zmax)

    # writes setup file
    lumi.setup(
        path=path,
        f=frequency,
        fs=frequencies,
        c0=c0,
        materials=materials,
        background=background,
        boundaries=boundaries,
        bbox=bbox,
        #
        frames=frames,
        sources=sources,
        monitors=monitors,
        #
        nres=nres,
        courant=courant,
        PML_decay=PML_decay,
        #
        views=views,
        shot_every=shot_every,
        tmax=tmax,
        shutoff=shutoff,
        device="gpu",
    )

setup(path, frequency, frequencies, w_line, h_metal, h_sub, l_stub, l_feed)
```

![diseño](img/stub_filter_layout.png)

```python
lumi.solve(path)
```

```python
sol = lumi.load(path)
m = sol["monitors"]
S21 = np.array(m["o2"]["modal"][0]["+"]) / np.array(m["o1"]["modal"][0]["-"])
i = np.argmin(abs(S21))
print(f"notch: {frequencies[i]:.2f} GHz, |S21| = {20 * np.log10(abs(S21[i])):.1f} dB")
```

```text
notch: 5.08 GHz, |S21| = -26.3 dB
```

