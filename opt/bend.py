import pylab
import EMpy
# import EMpy as em
import numpy
from functools import partial
from math import cos, pi, sin

import matplotlib.pyplot as plt
import numpy as np

from gdsfactory.cross_section import Section
import gdsfactory as gf
import bson
import utils
from utils import epscore, epsclad


def main(r, wwg, hwg, lwg=0, epscore=epscore, epsclad=epsclad):
    des = gf.Component()
    init = gf.Component()
    test = gf.Component()
    static = gf.Component()
    dev = gf.Component("bend")
    prob = dict()

    if lwg == 0:
        lwg = 2*wwg

    # Extrude the Path and the CrossSection
    wg1 = dev << gf.path.extrude(gf.path.straight(
        length=lwg), layer=(31, 0), width=wwg)
    wg2 = dev << gf.path.extrude(gf.path.straight(
        length=lwg), layer=(31, 0), width=wwg)
    stub = dev << gf.path.extrude(gf.path.straight(
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
    # des = dev << des

    wg1.connect("o2", des.ports["o1"])
    wg2.connect("o2", des.ports["o2"])
    stub.connect("o2", wg1.ports["o1"])

    dev.add_port("o1", port=wg1.ports["o1"])
    dev.add_port("o2", port=wg2.ports["o1"])

    static << des
    static << test
    static << dev
    # static.plot()
    # init.plot()
    utils.write_img("static", static, des_layer)
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
        for p in dev.ports.values()
    ]
    prob["ports"] = ports

    p = stub.ports["o1"]
    mode = utils.solve_modes(wwg, hwg, wm=.1, hm=.1,
                             dx=dx, epscore=epscore, epsclad=epsclad)[0]
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

    # dev.wr
    return dev


def ubend(r, wwg, hwg, l, lwg=0,
          λ=1.55, dx=0.05,
          hclad=0.1, hbase=0.1,
          epscore=epscore, epsclad=epsclad):
    des = gf.Component()
    init = gf.Component()
    test = gf.Component()
    static = gf.Component()
    dev = gf.Component("ubend")
    prob = dict()

    if lwg == 0:
        lwg = 2*wwg

    # Extrude the Path and the CrossSection
    wg1 = dev << gf.path.extrude(gf.path.straight(
        length=lwg), layer=(31, 0), width=wwg)
    wg2 = dev << gf.path.extrude(gf.path.straight(
        length=lwg), layer=(31, 0), width=wwg)

    lstub = lwg/2
    stub1 = dev << gf.path.extrude(gf.path.straight(
        length=lstub), layer=(51, 0), width=wwg)
    stub2 = dev << gf.path.extrude(gf.path.straight(
        length=lstub), layer=(51, 0), width=wwg)

    ld = l
    wd = 2*r+wwg
    des_layer = (100, 0)
    init_layer = (101, 0)
    des.add_polygon([(0, -wd/2), (ld, -wd/2), (ld, wd/2),
                    (0, wd/2)], layer=des_layer)
    init.add_polygon([(0, -wd/2), (ld, -wd/2), (ld, wd/2),
                     (0, wd/2)], layer=init_layer)
    des << init
    des.add_port(name="o1", center=(0, r), width=wwg,
                 orientation=180, layer=(100, 0))
    des.add_port(name="o2", center=(0, -r), width=wwg,
                 orientation=180, layer=(100, 0))
    # des = dev << des

    wg1.connect("o2", des.ports["o1"])
    wg2.connect("o2", des.ports["o2"])
    stub1.connect("o2", wg1.ports["o1"])
    stub2.connect("o2", wg2.ports["o1"])

    dev.add_port("o1", port=wg1.ports["o1"])
    dev.add_port("o2", port=wg2.ports["o1"])

    static << des
    static << test
    static << dev
    # static.plot()
    # init.plot()
    utils.write_img("static", static, hidden_layer=(des_layer, init_layer))
    utils.write_img("init", init, hidden_layer=(des_layer, ))

    # hbase = hclad = hwg/2
    ports = [
        {
            "center": p.center.flatten().tolist(),
            #  "lb": p.endpoints[0].flatten().tolist(),
            #  "ub": p.endpoints[1].flatten().tolist(),
            "normal": [cos(p.orientation/180*pi), sin(p.orientation/180*pi)],
            "endpoints": p.endpoints.tolist(),
        }
        for p in dev.ports.values()
    ]
    prob["ports"] = ports

    p = stub1.ports["o1"]
    modes, eps = utils.solve_modes(wwg, hwg, wm=.1, hm=.1,
                                   dx=dx, epscore=epscore, epsclad=epsclad)
    mode = modes[0]
    mode = {k: [np.real(mode[k]).tolist(), np.imag(mode[k]).tolist()]
            for k in mode}
    prob["signals"] = [
        {
            "wavelength": λ,
            "mode": mode,
            "center": p.center.flatten().tolist(),
            "width": eps.shape[0]*dx,
            "height": eps.shape[1]*dx,
            "eps": eps.tolist(),
            "normal": [cos(p.orientation/180*pi), sin(p.orientation/180*pi)],
            # "endpoints": p.endpoints.tolist(),
            # "ub": p.endpoints[1].flatten().tolist(),
            # "normal": (-p.normal).flatten().tolist()
        },
    ]

    prob["designs"] = [
        {
            "bbox": des.bbox.tolist(),
            # "bbox": des.bbox.flatten().tolist(),
        }
    ]

    prob["path_length"] = 2*(l+r+lwg+lstub)
    prob["hbase"] = hbase
    prob["hclad"] = hclad
    prob["hwg"] = hwg
    prob["epscore"] = epscore
    prob["epsclad"] = epsclad

    # bson.encode("prob.bson", prob)
    # np.savez
    # Encode the document to BSON
    bson_data = bson.dumps(prob)

    # Open a file for writing
    with open("prob.bson", "wb") as f:
        # Write the BSON data to the file
        f.write(bson_data)

    # dev.wr
    return prob


    # @save "$(@__DIR__)/layout.bson" static_mask signals ports designs dx λ ϵbase ϵclad ϵcore hsub hwg hclad modes l w
# if __module__
# main(1, .5, .22)
d = ubend(.4, .5, .22, 1.5)
1
