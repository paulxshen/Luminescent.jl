from copy import deepcopy
# from time import time
import datetime
import json
import subprocess
from functools import partial
from math import cos, pi, sin
import os
from .generic_tech import LAYER_STACK, LAYER, LAYER_VIEWS
import gdsfactory as gf
import bson
import scipy as sp
from .utils import *
from .layers import *
from .constants import *
from .utils import *
from .setup import *
from .sparams import *
from .materials import MATERIAL_LIBRARY


def inverse_design_problem(c,  lmin=.1, symmetries=[],
                           tparam_targets={}, sparam_targets={},
                           margin=XMARGIN,
                           maxiters=25, eta=.1, init=None,  # minloss=.01,
                           design_region_layer=LAYER.DESIGN,
                           #    design_guess_layer=LAYER.GUESS,
                           fill_layer=LAYER.WG,
                           void_layer=LAYER.WGCLAD,
                           layer_stack=LAYER_STACK,
                           plot=False, approx_2D=True, **kwargs):
    design_region_layer = tuple(design_region_layer)
    if not approx_2D:
        raise NotImplementedError(
            "3D inverse design feature must be requested from Luminescent AI info@luminescentai.com")

    if tparam_targets:
        target_type = "tparams"
        targets = tparam_targets
    elif sparam_targets:
        target_type = "sparams"
        targets = sparam_targets

    targets = {
        wl: {
            longname(k): v for k, v in d.items()
        } for wl, d in targets.items()}
    keys = sorted(set(sum([list(d.keys()) for d in targets.values()], [])))
    wavelengths = sorted(targets.keys())

    prob = sparams_problem(c, layer_stack=layer_stack, wavelengths=wavelengths,
                           study="inverse_design",
                           keys=keys, margin=margin, approx_2D=approx_2D, ** kwargs)
    prob["targets"] = targets
    prob["wavelengths"] = wavelengths
    prob["target_type"] = target_type
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

    prob["design_config"] = dict()
    l = get_layer(layer_stack, fill_layer)
    d = {"thickness": l.thickness, "material": l.material, "zmin": l.zmin}
    d["epsilon"] = MATERIAL_LIBRARY[d["material"]].epsilon
    d["ϵ"] = d["epsilon"]
    prob["design_config"]["fill"] = d

    l = get_layer(layer_stack, void_layer)
    d = {"thickness": l.thickness, "material": l.material, "zmin": l.zmin}
    d["epsilon"] = MATERIAL_LIBRARY[d["material"]].epsilon
    d["ϵ"] = d["epsilon"]
    prob["design_config"]["void"] = d

    prob["design_config"]["design_region_layer"] = design_region_layer
    prob["maxiters"] = maxiters
    return prob


def apply_design(c0,  sol):
    path = sol["path"]
    a = gf.Component()
    a.add_ref(c0)
    dl = sol["design_config"]["design_region_layer"]
    for i, d in enumerate(sol["designs"]):
        x0, y0 = d["bbox"][0]
        x1, y1 = d["bbox"][1]
        b = gf.Component()
        b.add_polygon(
            [(x0, y0), (x1, y0), (x1, y1), (x0, y1)], layer=dl)
        a = gf.boolean(a, b, "not", layer=dl)
    c = gf.Component()
    c << a
    for layer in c0.layers:
        if layer != dl:
            c.add_ref(c0.extract([layer]))
    # c.show()
    pic2gds(os.path.join(path, f"design{i+1}.png"), sol["dx"])
    g = gf.import_gds(os.path.join(path, f"design{i+1}.gds"), "TOP",
                      read_metadata=True)
    g = c << g
    g.drotate(90)
    g.xmin = x0
    g.ymin = y0
    return c
