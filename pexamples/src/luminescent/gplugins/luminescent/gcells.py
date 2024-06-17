from ast import mod
from multiprocessing.reduction import duplicate
import signal
import wave
from cv2 import add, norm
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
from .generic_tech import LAYER_STACK, LAYER_MAP, LAYER_VIEWS

import gdsfactory as gf
import bson
# import .utils as utils
from .utils import *
# from layers import *


def bend(r, wwg=.5, lwg=1, layer_map=LAYER_MAP, **kwargs):
    name = "bend"
    design = gf.Component()
    init = gf.Component()
    test = gf.Component()
    device = gf.Component()
    device = gf.Component("bend")
    prob = dict()

    # if lwg == 0:
    # lwg = 2*wwg
    lstub = lwg/2
    # Extrude the Path and the CrossSection
    wg1 = device << gf.path.extrude(gf.path.straight(
        length=lwg), layer=layer_map.WG, width=wwg)
    wg2 = device << gf.path.extrude(gf.path.straight(
        length=lwg), layer=layer_map.WG, width=wwg)
    arc = device << gf.path.extrude(gf.path.arc(
        r, start_angle=0), layer=layer_map.WG, width=wwg)

    ld = wd = 2*r
    design.add_polygon([(0, 0), (ld, 0), (ld, wd), (0, wd)],
                       layer=layer_map.DESIGN)
    design.add_port(name="o1", center=(0, r), width=wwg,
                    orientation=180, layer=layer_map.WG)
    design.add_port(name="o2", center=(r, 0), width=wwg,
                    orientation=-90, layer=layer_map.WG)
    design = device << design

    wg1.connect("o2", design.ports["o1"])
    wg2.connect("o2", design.ports["o2"])
    arc.connect("o1", wg2.ports["o2"])

    device.add_port("o1", port=wg1.ports["o1"])
    device.add_port("o2", port=wg2.ports["o1"])

    # utils.write_img("device", device, DESIGN)
    # utils.write_img("guess", init, DESIGN)
    c = add_bbox(device, layers=[layer_map.WGCLAD, layer_map.BOX],
                 nonport_margin=0.4)
    return c


# def ubend_template(r, l, wwg, lwg=0,
    #                wavelength=1.55, dx=0.05,
    #                #    ϵcore=ϵcore, ϵclad=ϵclad,
    #                dir="", init=None, lmin=0.1, layer_map=LAYER, **kwargs):
    # name = "ubend"
    # design = gf.Component("design")
    # device = gf.Component("ubend")
    # λ = wavelength

    # if lwg == 0:
    #     lwg = 4*wwg

    # # Extrude the Path and the CrossSection
    # wg1 = device << gf.path.extrude(gf.path.straight(
    #     length=lwg), layer=layer_map.WG, width=wwg)
    # wg2 = device << gf.path.extrude(gf.path.straight(
    #     length=lwg), layer=layer_map.WG, width=wwg)

    # ld = l
    # wd = 2*r+wwg

    # design.add_polygon([(0, -wd/2), (ld, -wd/2), (ld, wd/2),
    #                     (0, wd/2)], layer=DESIGN)
    # if init is not None:
    #     xmin = design.xmin
    #     ymin = design.ymin
    #     init = design << init
    #     init.xmin = xmin
    #     init.ymin = ymin
    # else:
    #     init = gf.Component()
    #     init.add_polygon([(0, -wd/2), (ld, -wd/2), (ld, wd/2),
    #                       (0, wd/2)], layer=GUESS)
    #     init = design << init

    # design.add_port(name="o1", center=(0, r), width=wwg,
    #                 orientation=180, layer=layer_map.PORT)
    # design.add_port(name="o2", center=(0, -r), width=wwg,
    #                 orientation=180, layer=layer_map.PORT)
    # # design = device << design

    # wg1.connect("o2", design.ports["o1"], allow_layer_mismatch=True)
    # wg2.connect("o2", design.ports["o2"], allow_layer_mismatch=True)
    # # stub1.connect("o2", wg1.ports["o1"])
    # # stub2.connect("o2", wg2.ports["o1"])

    # device.add_port("o1", port=wg1.ports["o1"])
    # device.add_port("o2", port=wg2.ports["o1"])

    # # design = device << design
    # device.add_ref(design, name="design")
    # # device << gf.components.re
    # # device << test
    # device.show()
    # # init.show()
    # sources = {
    #     "o1": {
    #         "port": "o1",
    #         "wavelengths": [λ],
    #         "mode_numbers": [0],
    #         "powers": [1],
    #         "width": wwg,
    #         "endpoints": device.ports["o1"].endpoints.tolist(),
    #         "center": (device.ports["o1"].center),
    #         "orientation": device.ports["o1"].orientation,
    #     }}
    # for k in sources:
    #     a = sources[k]["orientation"]/180*pi
    #     normal = [cos(a), sin(a)]
    #     sources[k]["normal"] = normal
    #     sources[k]["center"] -= .001*np.array(normal)
    #     sources[k]["center"] = sources[k]["center"].flatten().tolist()

    # monitors = {
    #     "o2":        {
    #         "powers": [1],
    #         "wavelengths": [λ],
    #         "mode_numbers": [0],
    #         "port": "o2"},
    #     # {"power": 1,"wavelength": λ,"port": "o1"},
    # }
    # # device.metadata["sources"] = sources
    # # device.metadata["monitors"] = monitors
    # return device, sources, monitors, design
