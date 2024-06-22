from copy import deepcopy
# from time import time
import datetime
import json
import subprocess
from functools import partial
from math import cos, pi, sin
import os
import matplotlib.pyplot as plt
import numpy as np

from gdsfactory.cross_section import Section
from sortedcontainers import SortedDict, SortedSet
from .generic_tech import LAYER_STACK, LAYER_MAP, LAYER_VIEWS
import gdsfactory as gf
import bson
import scipy as sp
from .utils import *
from .layers import *
from .constants import *
from .utils import *
from .setup import *
from .sparams import *


def inverse_design_problem(c, sparam_targets, lmin=.1, symmetries=[],
                           maxiters=25,
                           design_region_layer=LAYER_MAP.DESIGN,
                           design_guess_layer=LAYER_MAP.GUESS,
                           design_layer=LAYER_MAP.WG,
                           void_layer=LAYER_MAP.WGCLAD,
                           layer_stack=LAYER_STACK,
                           plot=False, **kwargs):
    targets = {
        wl: {longname(k): v for k, v in d.items()}for wl, d in sparam_targets.items()}
    keys = list(targets.values())[0].keys()
    prob = sparams_problem(c, layer_stack=layer_stack, keys=keys, **kwargs)
    prob["study"] = "inverse_design"
    polys = c.extract([design_region_layer]).get_polygons()
    if not symmetries:
        symmetries = [[]]*len(polys)
    else:
        if type(symmetries[0]) is int:
            symmetries = [symmetries]*len(polys)

    def _bbox(b):
        return [[b.left/1e3, b.bottom/1e3], [b.right/1e3, b.top/1e3]]
    prob["designs"] = [
        {
            "layer": design_layer,
            "bbox": _bbox(p.bbox()),
            "symmetries": s,
            "guess": None,
            "lmin": lmin,
        } for p, s in zip(list(polys.values())[0], symmetries)
    ]

    prob["design_config"] = dict()
    # prob[""]
    l = get_layer(layer_stack, design_layer)
    d = {"thickness": l.thickness, "material": l.material, "zmin": l.zmin}
    # d = vars(prob["design_layer"])
    d["ϵ"] = MATERIAL_EPS[d["material"]]
    prob["design_config"]["fill"] = d

    l = get_layer(layer_stack, void_layer)
    d = {"thickness": l.thickness, "material": l.material, "zmin": l.zmin}
    d["ϵ"] = MATERIAL_EPS[d["material"]]
    prob["design_config"]["void"] = d

    prob["design_layer"] = d
    prob["targets"] = targets
    prob["maxiters"] = maxiters
    return prob


def apply_design(c0, prob, sol):
    path = prob["path"]
    a = gf.Component()
    a.add_ref(c0)
    dl = 0
    for i, d in enumerate(prob["designs"]):
        dl = d["layer"]
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
    c.show()
    pic2gds(os.path.join(path, f"design{i+1}.png"), prob["dx"])
    g = gf.import_gds(os.path.join(path, f"design{i+1}.gds"), "TOP",
                      read_metadata=True)
    g = c << g
    g.drotate(90)
    g.xmin = x0
    g.ymin = y0
    return c
