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

pdk = get_generic_pdk()
pdk.activate()
LAYER_VIEWS = pdk.layer_views
# include_layers=[layer_map.WG, layer_map.WAFER], ** kwargs):


def setup(c, study, name, dx, wavelengths=[], signals=[], layer_stack=LAYER_STACK, exclude_layers=[DESIGN, GUESS], ):
    prob = dict()
    prob["name"] = name
    ports = {
        k: {
            "center": p.center.flatten().tolist(),
            #  "lb": p.endpoints[0].flatten().tolist(),
            #  "ub": p.endpoints[1].flatten().tolist(),
            "normal": [cos(p.orientation/180*pi), sin(p.orientation/180*pi)],
            "endpoints": p.endpoints.tolist(),
            # **targets[k],
        }
        for k, p in zip(c.ports, c.ports.values(), ) if k.startswith("o")
    }
    prob["ports"] = ports

    mode_solutions = []
    if study == "sparams":
        1
    else:
        1
        # prob["signals"] = c.metadata["signals"]
        # prob["targets"] = c.metadata["targets"]

    hcore = layer_stack.layers["core"].thickness
    zcenter = layer_stack.layers["core"].zmin+hcore/2
    for s in signals.values():
        for wl in s["wavelengths"]:
            duplicate = False
            wwg = s["width"]
            for m in mode_solutions:
                if m["wavelength"] == wl and m["width"] == s["width"]:
                    m["ports"].append(s["port"])
                    duplicate = True
            if not duplicate:
                wm = .2
                hm = .2
                w = wwg+2*wm
                h = hcore+2*hm
                a = []
                layers = set(c.layers)-set(exclude_layers)
                layers = sorted(layers, key=lambda layer: next(
                    x for x in layer_stack.layers.values() if x.layer == layer).mesh_order)
                a = None
                epsmin = 100

                layer_views = LAYER_VIEWS.copy()
                layer_views.layer_views["WGCLAD"].visible = True

                for layer in layers:
                    layer_views.layer_views[f"{layer}"] = LAYER_VIEWS.layer_views["WGCLAD"]
                    layer_views.layer_views[f"{layer}"].layer = layer

                    # for layer in [layer_map.WAFER]:
                    print(layer)
                    eps = material_eps[next(
                        x for x in layer_stack.layers.values() if x.layer == layer).material]
                    mask = raster_slice(c.to_3d(
                        layer_stack=layer_stack,
                        layer_views=layer_views,
                        exclude_layers=list(set(c.layers)-set([layer]))), dx, w=w, h=h, center=s["center"]+[zcenter], normal=s["normal"]+[0])
                    if eps < epsmin:
                        epsmin = eps
                    _a = eps * mask
                    if mask is not None:
                        if a is None:
                            a = _a.copy()
                        else:
                            # a += _a*np.where(a == 0, 1, 0)
                            a += _a*np.where(a == 0, 1, 0)
                eps = a+np.where(a == 0, 1, 0)*epsmin
                plt.clf()
                plt.imshow(eps)
                plt.show()
                _modes = utils.solve_modes(eps, λ=wl, dx=dx,)
                # mode = _modes[0]
                modes = [{k: [np.real(mode[k]).tolist(), np.imag(
                    mode[k]).tolist()] for k in mode} for mode in _modes]
                mode_solutions.append({
                    "modes": modes,
                    "wavelength": wl,
                    "port": s["port"],
                    "size": [w, h],
                    "eps": eps.tolist(),
                })
    prob["mode_solutions"] = mode_solutions

    device_svg = utils.write_img("device", c,
                                 hidden_layer=(DESIGN, GUESS))
    prob["masks"] = {
        "device": {
            "svg": device_svg,
            "bbox": c.bbox.tolist(),
        }}
    prob["dx"] = dx
    λc = (wavelengths[0]+wavelengths[-1])/2
    prob["λc"] = λc
    prob["wavelengths"] = wavelengths
    prob["mode_solutions"] = mode_solutions
    # prob["path_length"] = 2*(l+r+lwg)

    m = round(λc/2/dx)*dx
    prob["margins"] = [[m, m], [m, m]]
    return prob


def inverse_design_problem(c,  lmin, dx, signals, targets, design, name="", layer_stack=LAYER_STACK, **kwargs):
    # wavelengths = c.metadata["signals"]["o1"]["wavelengths"]
    wavelengths = signals["o1"]["wavelengths"]
    prob = setup(c, name=name, study="inverse_design",
                 dx=dx, wavelengths=wavelengths, signals=signals, layer_stack=layer_stack)
    prob["lmin"] = lmin

    prob["targets"] = targets
    prob["signals"] = signals
    designs = [
        {
            # "bbox": c.named_references["design"].bbox(),
            "bbox": design.bbox.flatten().tolist(),
            "symmetry_dims": [2],
        }
    ]
    prob["designs"] = designs

    # design = c.child["design"]
    init_svg = utils.write_img(
        "guess", design, hidden_layer=(DESIGN, ))
    prob["masks"]["guess"] = {
        "svg": init_svg,
        "bbox": design.bbox.tolist(),
    }

    # Encode the document to BSON

    return prob


def solve(prob,):
    bson_data = bson.dumps(prob)

    # # Open a file for writing
    with open(os.path.join("_prob.bson"), "wb") as f:
        # Write the BSON data to the file
        f.write(bson_data)

    # @save "$(@__DIR__)/layout.bson" device_mask signals ports designs dx λ ϵbase ϵclad ϵcore hbase hwg hclad modes l w
# if __module__
# main(1, .5, .22)
