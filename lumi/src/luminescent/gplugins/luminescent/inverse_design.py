from .materials import MATERIALS
from .sparams import *
from .setup import *
from .constants import *
from .layers import *
from .utils import *
import scipy as sp
import bson
import gdsfactory as gf
from copy import deepcopy
# from time import time
import datetime
import json
import subprocess
from functools import partial
from math import cos, pi, sin
import os
from gdsfactory.generic_tech import LAYER_STACK, LAYER


def inverse_design_problem(c,  targets=dict(), preset=None,
                           lmin=.1, symmetries=[],
                           weights=dict(),
                           iters=25, eta=2., init=None,  # minloss=.01,
                           design_region_layer=DESIGN_LAYER,
                           #    design_guess_layer=LAYER.GUESS,
                           fill_layer=LAYER.WG,
                           void_layer=None,
                           layer_stack=LAYER_STACK, materials=MATERIALS,
                           contrast=20,
                           plot=False, approx_2D=True,
                           restart=True, lowmem=False, **kwargs):
    design_region_layer = tuple(design_region_layer)
    # if not approx_2D:
    #     raise NotImplementedError(
    #         "3D inverse design feature must be requested from Luminescent AI info@luminescentai.com")

    if preset is not None:
        if preset["name"] == "phase_shifter":
            keys = ["o2@0,o1@0"]
            wavelengths = [preset["wavelength"]]
    else:
        targets1 = {}
        keys = SortedSet()
        wavelengths = SortedSet()
        for k, d in targets.items():
            d = {
                wl: {
                    longname(k): v for k, v in d.items()
                } for wl, d in d.items()}
            targets1[k] = d
            if k in ["sparams", "tparams"]:
                for s in sum([list(d.keys()) for d in d.values()], []):
                    keys.add(s)
            for wl in d:
                wavelengths.add(wl)
        targets = targets1
        keys = list(keys)
        wavelengths = list(wavelengths)

    prob = sparams_problem(c,
                           layer_stack=layer_stack, materials=materials,
                           wavelengths=wavelengths,
                           study="inverse_design",
                           keys=keys,
                           approx_2D=approx_2D, ** kwargs)
    prob["restart"] = restart
    prob["weights"] = {**{
        "tparams": 1,
        "sparams": 1,
        "phasediff": 1,
    }, **weights}
    prob["lowmem"] = lowmem
    prob["preset"] = preset
    prob["targets"] = targets
    prob["wavelengths"] = wavelengths
    prob["contrast"] = contrast
    # prob["init"] = init
    prob = {**prob, **kwargs}
    prob["eta"] = eta
    polys = c.extract([design_region_layer]).get_polygons()
    if not symmetries:
        symmetries = [[]]*len(polys)
    else:
        if type(symmetries[0]) is not list:
            symmetries = [symmetries]*len(polys)

    def _bbox(b):
        return [[b.left/1e3, b.bottom/1e3], [b.right/1e3, b.top/1e3]]
    prob["designs"] = [
        {
            "layer": design_region_layer,
            "bbox": _bbox(p.bbox()),
            "symmetries": s,
            "init": init,
            "lmin": lmin,
        } for p, s in zip(list(polys.values())[0], symmetries)
    ]
    epsmin = np.min(prob["eps_2D"])

    prob["design_config"] = dict()
    l = get_layers(layer_stack, fill_layer)[0]
    d = {"thickness": l.thickness, "material": l.material, "zmin": l.zmin}
    d["epsilon"] = materials[d["material"]].epsilon
    d["ϵ"] = d["epsilon"]
    d["layer"] = fill_layer
    prob["design_config"]["fill"] = d

    if void_layer is not None:
        l = get_layers(layer_stack, void_layer)[0]
        d = {"thickness": l.thickness, "material": l.material, "zmin": l.zmin}
        d["layer"] = void_layer
        d["epsilon"] = materials[d["material"]].epsilon
        d["ϵ"] = d["epsilon"]
    else:
        d = copy.deepcopy(d)
        d["epsilon"] = epsmin
        d["ϵ"] = d["epsilon"]
    prob["design_config"]["void"] = d

    prob["design_config"]["design_region_layer"] = design_region_layer
    prob["iters"] = iters
    return prob


def apply_design(c0,  sol):
    path = sol["path"]
    a = gf.Component()
    a.add_ref(c0)
    fill = sol["design_config"]["fill"]["layer"]
    dl = sol["design_config"]["design_region_layer"]
    for i, d in enumerate(sol["designs"]):
        x0, y0 = d["bbox"][0]
        x1, y1 = d["bbox"][1]
        b = gf.Component()
        b.add_polygon(
            [(x0, y0), (x1, y0), (x1, y1), (x0, y1)], layer=fill)
        a = gf.boolean(a, b, "not", layer=fill)
    c = gf.Component()
    c << a
    # for layer in c0.layers:
    #     if layer != dl:
    #         c.add_ref(c0.extract([layer]))
    # c.show()
    # raise ValueError()
    pic2gds(os.path.join(path, f"design{i+1}.png"), sol["dx"])
    g = gf.import_gds(os.path.join(path, f"design{i+1}.gds"), "TOP",
                      read_metadata=True)
    polygons = g.get_polygons(merge=True)
    g = gf.Component()
    for p in polygons[1]:
        g.add_polygon(p, layer=fill)
    g = c << g
    g.drotate(90)
    g.xmin = x0
    g.ymin = y0
    return c
