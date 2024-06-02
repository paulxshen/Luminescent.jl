from ast import mod
from multiprocessing.reduction import duplicate
import signal
import wave
from cv2 import norm
import pylab
import EMpy
# import EMpy as em
import numpy
from functools import partial
from math import cos, pi, sin
import os
import matplotlib.pyplot as plt
import numpy as np

from gdsfactory.cross_section import Section
from gdsfactory.generic_tech import LAYER_STACK, get_generic_pdk, LAYER
import gdsfactory as gf
import bson
import utils
from utils import *
from layers import *


# # def bend(r, wwg, hwg, lwg=0, ϵcore=ϵcore, ϵclad=ϵclad, lmin=0.1, **kwargs):
#     name = "bend"
#     design = gf.Component()
#     init = gf.Component()
#     test = gf.Component()
#     device = gf.Component()
#     device = gf.Component("bend")
#     prob = dict()

#     if lwg == 0:
#         lwg = 2*wwg
#     lstub = lwg/2
#     # Extrude the Path and the CrossSection
#     wg1 = device << gf.path.extrude(gf.path.straight(
#         length=lwg), layer=(31, 0), width=wwg)
#     wg2 = device << gf.path.extrude(gf.path.straight(
#         length=lwg), layer=(31, 0), width=wwg)
#     stub = device << gf.path.extrude(gf.path.straight(
#         length=lstub), layer=(51, 0), width=wwg)
#     arc = init << gf.path.extrude(gf.path.arc(
#         r, start_angle=0), layer=(101, 0), width=wwg)

#     ld = wd = 2*r
#     DESIGN = (100, 0)
#     design.add_polygon([(0, 0), (ld, 0), (ld, wd), (0, wd)], layer=DESIGN)
#     design.add_port(name="o1", center=(0, r), width=wwg,
#                     orientation=180, layer=(100, 0))
#     design.add_port(name="o2", center=(r, 0), width=wwg,
#                     orientation=-90, layer=(100, 0))
#     # design = device << design

#     wg1.connect("o2", design.ports["o1"])
#     wg2.connect("o2", design.ports["o2"])
#     stub.connect("o2", wg1.ports["o1"])

#     device.add_port("o1", port=wg1.ports["o1"])
#     device.add_port("o2", port=wg2.ports["o1"])

#     device << design
#     device << test
#     device << device
#     # device.show()
#     # init.show()
#     utils.write_img("device", device, DESIGN)
#     utils.write_img("guess", init, DESIGN)

#     dx = 0.05
#     λ = 1.55
#     # hbase = hclad = hwg/2
#     ports = [
#         {
#             "center": p.center.flatten().tolist(),
#             #  "lb": p.endpoints[0].flatten().tolist(),
#             #  "ub": p.endpoints[1].flatten().tolist(),
#             "normal": [cos(p.orientation/180*pi), sin(p.orientation/180*pi)],
#             "endpoints": p.endpoints.tolist(),
#         }
#         for p in device.ports.values()
#     ]
#     prob["ports"] = ports

#     p = stub.ports["o1"]
#     modes, ϵ = utils.solve_modes(wwg, hwg, wm=.1, hm=.1,
#                                  dx=dx, ϵcore=ϵcore, ϵclad=ϵclad)[0]
#     mode = modes[0]
#     mode = {k: [np.real(mode[k]).tolist(), np.imag(mode[k]).tolist()]
#             for k in mode}
#     prob["sources"] = [
#         {
#             "wavelength": λ,
#             "mode": mode,
#             "center": p.center.flatten().tolist(),
#             "endpoints": p.endpoints.tolist(),
#             # "ub": p.endpoints[1].flatten().tolist(),
#             "normal": (-p.normal).flatten().tolist()},
#     ]

#     prob["design"] = [
#         {
#             "bbox": design.bbox.flatten().tolist(),
#             # "bbox": design.bbox.flatten().tolist(),
#         }
#     ]
#     # bson.encode("prob.bson", prob)
#     # np.savez
#     # Encode the document to BSON
#     bson_data = bson.dumps(prob)

#     # Open a file for writing
#     with open("prob.bson", "wb") as f:
#         # Write the BSON data to the file
#         f.write(bson_data)

#     # device.wr
#     return device


def ubend_template(r, l, wwg, lwg=0,
                   wavelength=1.55, dx=0.05,
                   #    ϵcore=ϵcore, ϵclad=ϵclad,
                   dir="", init=None, lmin=0.1, layer_map=LAYER, **kwargs):
    name = "ubend"
    design = gf.Component("design")
    device = gf.Component("ubend")
    λ = wavelength

    if lwg == 0:
        lwg = 4*wwg

    # Extrude the Path and the CrossSection
    wg1 = device << gf.path.extrude(gf.path.straight(
        length=lwg), layer=layer_map.WG, width=wwg)
    wg2 = device << gf.path.extrude(gf.path.straight(
        length=lwg), layer=layer_map.WG, width=wwg)

    ld = l
    wd = 2*r+wwg

    design.add_polygon([(0, -wd/2), (ld, -wd/2), (ld, wd/2),
                        (0, wd/2)], layer=DESIGN)
    if init is not None:
        xmin = design.xmin
        ymin = design.ymin
        init = design << init
        init.xmin = xmin
        init.ymin = ymin
    else:
        init = gf.Component()
        init.add_polygon([(0, -wd/2), (ld, -wd/2), (ld, wd/2),
                          (0, wd/2)], layer=GUESS)
        init = design << init

    design.add_port(name="o1", center=(0, r), width=wwg,
                    orientation=180, layer=layer_map.PORT)
    design.add_port(name="o2", center=(0, -r), width=wwg,
                    orientation=180, layer=layer_map.PORT)
    # design = device << design

    wg1.connect("o2", design.ports["o1"], allow_layer_mismatch=True)
    wg2.connect("o2", design.ports["o2"], allow_layer_mismatch=True)
    # stub1.connect("o2", wg1.ports["o1"])
    # stub2.connect("o2", wg2.ports["o1"])

    device.add_port("o1", port=wg1.ports["o1"])
    device.add_port("o2", port=wg2.ports["o1"])

    # design = device << design
    device.add_ref(design, name="design")
    # device << gf.components.re
    # device << test
    device.show()
    # init.show()
    sources = {
        "o1": {
            "port": "o1",
            "wavelengths": [λ],
            "mode_numbers": [0],
            "powers": [1],
            "width": wwg,
            "endpoints": device.ports["o1"].endpoints.tolist(),
            "center": (device.ports["o1"].center),
            "orientation": device.ports["o1"].orientation,
        }}
    for k in sources:
        a = sources[k]["orientation"]/180*pi
        normal = [cos(a), sin(a)]
        sources[k]["normal"] = normal
        sources[k]["center"] -= .001*np.array(normal)
        sources[k]["center"] = sources[k]["center"].flatten().tolist()

    targets = {
        "o2":        {
            "powers": [1],
            "wavelengths": [λ],
            "mode_numbers": [0],
            "port": "o2"},
        # {"power": 1,"wavelength": λ,"port": "o1"},
    }
    # device.metadata["sources"] = sources
    # device.metadata["targets"] = targets
    return device, sources, targets, design
