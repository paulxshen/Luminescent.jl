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
import gdsfactory as gf
import bson
import utils
from utils import *
from layers import *


def bend(r, wwg, hwg, lwg=0, ϵcore=ϵcore, ϵclad=ϵclad, lmin=0.1, **kwargs):
    name = "bend"
    des = gf.Component()
    init = gf.Component()
    test = gf.Component()
    device = gf.Component()
    device = gf.Component("bend")
    prob = dict()

    if lwg == 0:
        lwg = 2*wwg
    lstub = lwg/2
    # Extrude the Path and the CrossSection
    wg1 = device << gf.path.extrude(gf.path.straight(
        length=lwg), layer=(31, 0), width=wwg)
    wg2 = device << gf.path.extrude(gf.path.straight(
        length=lwg), layer=(31, 0), width=wwg)
    stub = device << gf.path.extrude(gf.path.straight(
        length=lstub), layer=(51, 0), width=wwg)
    arc = init << gf.path.extrude(gf.path.arc(
        r, start_angle=0), layer=(101, 0), width=wwg)

    ld = wd = 2*r
    des_layer = (100, 0)
    des.add_polygon([(0, 0), (ld, 0), (ld, wd), (0, wd)], layer=des_layer)
    des.add_port(name="o1", center=(0, r), width=wwg,
                 orientation=180, layer=(100, 0))
    des.add_port(name="o2", center=(r, 0), width=wwg,
                 orientation=-90, layer=(100, 0))
    # des = device << des

    wg1.connect("o2", des.ports["o1"])
    wg2.connect("o2", des.ports["o2"])
    stub.connect("o2", wg1.ports["o1"])

    device.add_port("o1", port=wg1.ports["o1"])
    device.add_port("o2", port=wg2.ports["o1"])

    device << des
    device << test
    device << device
    # device.show()
    # init.show()
    utils.write_img("device", device, des_layer)
    utils.write_img("init", init, des_layer)

    dx = 0.05
    λ = 1.55
    # hbase = hclad = hwg/2
    ports = [
        {
            "center": p.center.flatten().tolist(),
            #  "lb": p.endpoints[0].flatten().tolist(),
            #  "ub": p.endpoints[1].flatten().tolist(),
            "normal": [cos(p.orientation/180*pi), sin(p.orientation/180*pi)],
            "endpoints": p.endpoints.tolist(),
        }
        for p in device.ports.values()
    ]
    prob["ports"] = ports

    p = stub.ports["o1"]
    modes, ϵ = utils.solve_modes(wwg, hwg, wm=.1, hm=.1,
                                 dx=dx, ϵcore=ϵcore, ϵclad=ϵclad)[0]
    mode = modes[0]
    mode = {k: [np.real(mode[k]).tolist(), np.imag(mode[k]).tolist()]
            for k in mode}
    prob["signals"] = [
        {
            "wavelength": λ,
            "mode": mode,
            "center": p.center.flatten().tolist(),
            "endpoints": p.endpoints.tolist(),
            # "ub": p.endpoints[1].flatten().tolist(),
            "normal": (-p.normal).flatten().tolist()},
    ]

    prob["des"] = [
        {
            "bbox": des.bbox.flatten().tolist(),
            # "bbox": des.bbox.flatten().tolist(),
        }
    ]
    # bson.encode("prob.bson", prob)
    # np.savez
    # Encode the document to BSON
    bson_data = bson.dumps(prob)

    # Open a file for writing
    with open("prob.bson", "wb") as f:
        # Write the BSON data to the file
        f.write(bson_data)

    # device.wr
    return device


def ubend(r, wwg, hwg, l, lwg=0,
          λ=1.55, dx=0.05,
          hclad=0.1, hbase=0.1,
          ϵcore=ϵcore, ϵclad=ϵclad,
          dir="", init=None, lmin=0.1, **kwargs):
    name = "ubend"
    des = gf.Component()
    device = gf.Component("ubend")
    prob = dict()

    if lwg == 0:
        lwg = 4*wwg

    # Extrude the Path and the CrossSection
    wg1 = device << gf.path.extrude(gf.path.straight(
        length=lwg), layer=(31, 0), width=wwg)
    wg2 = device << gf.path.extrude(gf.path.straight(
        length=lwg), layer=(31, 0), width=wwg)

    # lstub = lwg/2
    # stub1 = device << gf.path.extrude(gf.path.straight(
    #     length=lstub), layer=(51, 0), width=wwg)
    # stub2 = device << gf.path.extrude(gf.path.straight(
    #     length=lstub), layer=(51, 0), width=wwg)

    ld = l
    wd = 2*r+wwg

    des.add_polygon([(0, -wd/2), (ld, -wd/2), (ld, wd/2),
                    (0, wd/2)], layer=des_layer)
    if init is not None:
        xmin = des.xmin
        ymin = des.ymin
        init = des << init
        init.xmin = xmin
        init.ymin = ymin
    else:
        init = gf.Component()
        init.add_polygon([(0, -wd/2), (ld, -wd/2), (ld, wd/2),
                          (0, wd/2)], layer=init_layer)
        init = des << init

    des.add_port(name="o1", center=(0, r), width=wwg,
                 orientation=180, layer=(100, 0))
    des.add_port(name="o2", center=(0, -r), width=wwg,
                 orientation=180, layer=(100, 0))
    # des = device << des

    wg1.connect("o2", des.ports["o1"])
    wg2.connect("o2", des.ports["o2"])
    # stub1.connect("o2", wg1.ports["o1"])
    # stub2.connect("o2", wg2.ports["o1"])

    device.add_port("o1", port=wg1.ports["o1"])
    device.add_port("o2", port=wg2.ports["o1"])

    device << des
    # device << test
    device.show()
    # init.show()

    device_svg = utils.write_img(dir, "device", device,
                                 hidden_layer=(des_layer, init_layer))
    init_svg = utils.write_img(dir, "init", des, hidden_layer=(des_layer, ))

    # hbase = hclad = hwg/2
    prob["name"] = "ubend"
    targets = {
        "o1":        {"power": 1},
        "o2":   {"power": 1},
    }
    ports = {
        k: {
            "center": p.center.flatten().tolist(),
            #  "lb": p.endpoints[0].flatten().tolist(),
            #  "ub": p.endpoints[1].flatten().tolist(),
            "normal": [cos(p.orientation/180*pi), sin(p.orientation/180*pi)],
            "endpoints": p.endpoints.tolist(),
            **targets[k],
        }
        for k, p in zip(device.ports, device.ports.values(), )
    }
    prob["ports"] = ports

    p = wg1.ports["o1"]
    modes, ϵ, L = utils.solve_modes(wwg, hwg, wm=.2, hm=.2,
                                    dx=dx, ϵcore=ϵcore, ϵclad=ϵclad)
    mode = modes[0]
    mode = {k: [np.real(mode[k]).tolist(), np.imag(mode[k]).tolist()]
            for k in mode}
    width = L[0]
    height = L[1]
    prob["signals"] = {
        "o1": {
            "wavelength": λ,
            "mode": mode,
            "center": p.center.flatten().tolist(),
            "width": width,
            "height": height,
            "size": [width, height],
            "endpoints": [
                ((p.endpoints[0]-p.center) /
                 np.linalg.norm(p.endpoints[0]-p.center)*width/2).tolist(),
                ((p.endpoints[1]-p.center) /
                 np.linalg.norm(p.endpoints[1]-p.center)*width/2).tolist(),
            ],
            "wwg": wwg,
            "ϵ": ϵ.tolist(),
            "normal": [-cos(p.orientation/180*pi), -sin(p.orientation/180*pi)],
            # "endpoints": p.endpoints.tolist(),
            # "ub": p.endpoints[1].flatten().tolist(),
            # "normal": (-p.normal).flatten().tolist()
        },
    }

    prob["designs"] = [
        {
            "bbox": des.bbox.tolist(),
            # "bbox": des.bbox.flatten().tolist(),
            "symmetry_dims": [2],
        }
    ]

    prob["layers"] = {
        "device": {
            "svg": device_svg,
            "bbox": device.bbox.tolist(),
        },
        "init": {
            "svg": init_svg,
            "bbox": init.bbox.tolist(),
        }
    }
    prob["dx"] = dx
    prob["λ"] = λ
    prob["path_length"] = 2*(l+r+lwg)
    prob["hbase"] = hbase
    prob["hclad"] = hclad
    prob["hwg"] = hwg
    prob["ϵcore"] = ϵcore
    prob["ϵclad"] = ϵclad
    prob["ϵbase"] = ϵclad
    prob["targets"] = targets
    prob["lmin"] = lmin

    m = round( λ/2/dx)*dx
    prob["margins"] = [[0, m], [m, m]]
    # bson.encode("prob.bson", prob)
    # np.savez
    # Encode the document to BSON
    bson_data = bson.dumps(prob)

    # Open a file for writing
    with open(os.path.join(dir, "prob.bson"), "wb") as f:
        # Write the BSON data to the file
        f.write(bson_data)

    # device.wr
    # return prob
    finish(device, name)
    return device

    # @save "$(@__DIR__)/layout.bson" device_mask signals ports designs dx λ ϵbase ϵclad ϵcore hbase hwg hclad modes l w
# if __module__
# main(1, .5, .22)
