from copy import deepcopy
# from time import time
import datetime
import json
import subprocess
import sys
import time
from functools import partial
from math import cos, pi, sin
import os
import numpy as np

from gdsfactory.cross_section import Section
from sortedcontainers import SortedDict, SortedSet
from sympy import N
from .generic_tech import LAYER_STACK, LAYER, LAYER_VIEWS
import gdsfactory as gf
import bson
from .utils import *
from .layers import *
from .constants import *
from .utils import *

# include_layers=[LAYER.WG, LAYER.WAFER], ** kwargs):


def setup(c, study, center_wavelength,  dx,
          zmargin,
          xmargin, name="",
          runs=[], wavelengths=[], sources=[], layer_stack=LAYER_STACK, layer_views=LAYER_VIEWS, exclude_layers=[
              LAYER.DESIGN, GUESS], approx_2D=False,
          gpu=None,
          path=PATH, plot=False, **kwargs):
    prob = dict()
    prob["path"] = os.path.join(
        path, datetime.datetime.now().isoformat(timespec="seconds").replace(":", "-"))
    prob["name"] = name
    prob["gpu_backend"] = gpu if gpu else ""
    ports = {
        p.name: {
            "center": p.center,
            #  "lb": p.endpoints[0].flatten().tolist(),
            #  "ub": p.endpoints[1].flatten().tolist(),
            "normal": [cos(p.orientation/180*pi), sin(p.orientation/180*pi)],
            # "endpoints": p.endpoints.tolist(),
            # **monitors[k],
        }
        for p in c.get_ports_list(prefix="o")
    }
    prob["ports"] = ports

    mode_solutions = []
    if study == "detailed_sparams":
        1
    else:
        1
        # prob["sources"] = c.metadata["sources"]
        # prob["monitors"] = c.metadata["monitors"]
    layers = set(c.layers)-set(exclude_layers)
    layers = sorted(layers, key=lambda layer: -get_layer(
        layer_stack, layer).mesh_order)

    hcore = layer_stack.layers["core"].thickness
    zcenter = layer_stack.layers["core"].zmin+hcore/2
    h = hcore+2*zmargin
    prob["mode_height"] = h
    l, w = c.bbox_np()[1]-c.bbox_np()[0]
    center = [c.bbox_np()[0][0], c.bbox_np()[0][1]+w/2, zcenter]
    normal = [1, 0, 0]
    eps = material_voxelate(c, dx, center, l, w, h,
                            normal, layers, layer_stack)
    eps_2D = eps[:, :, int(eps.shape[2]/2)]
    prob["eps_3D"] = eps.tolist()
    prob["eps_2D"] = eps_2D.tolist()
    prob["zmin"] = -zmargin
    prob["d"] = 2 if approx_2D else 3

    for run in runs:
        for s in list(run["sources"].values())+list(run["monitors"].values()):
            for wl in s["wavelength_mode_numbers"]:
                duplicate = False
                w = s["width"]
                mode_numbers = s["wavelength_mode_numbers"][wl]
                for m in mode_solutions:
                    if m["wavelength"] == wl and max(mode_numbers) < len(m["modes"]) and m["width"] == s["width"]:
                        m["ports"].append(s["port"])
                        duplicate = True
                if not duplicate:
                    normal = s["normal"]+[0]
                    center = list(s["center"])+[zcenter]

                    center = np.array(center)-.001*np.array(normal)
                    eps = material_slice(
                        c, dx, center, w, h, normal, layers, layer_stack, layer_views=layer_views)
                    _modes = solve_modes(
                        eps, λ=wl, dx=dx, neigs=max(mode_numbers)+1, plot=plot)
                    # _modes = [0]
                    # mode = _modes[0]
                    modes = [{k: [np.real(mode[k]).tolist(), np.imag(
                        mode[k]).tolist()] for k in mode} for mode in _modes]
                    # for mn in mode_numbers:
                    mode_solutions.append({
                        "modes": modes,
                        "wavelength": wl,
                        # "mode_number": mn,
                        "ports": [s["port"]],
                        "size": [w, h],
                        "eps": eps.tolist(),
                        "width": w,
                        "zcenter": zcenter,
                    })
    prob["mode_solutions"] = mode_solutions
    prob["runs"] = runs
    # device_svg = write_img("device", c,
    #                        hidden_layer=(DESIGN, GUESS))
    prob["components"] = {
        "device": {
            # "svg": device_svg,
            "bbox": c.bbox_np().tolist(),
        }}
    prob["dx"] = dx
    prob["component"] = c
    # λc = (wavelengths[0]+wavelengths[-1])/2
    if center_wavelength is None:
        center_wavelength = mode_solutions[0]["wavelength"]
    prob["λc"] = center_wavelength
    # prob["wavelengths"] = wavelengths
    prob["mode_solutions"] = mode_solutions
    # prob["path_length"] = 2*(l+r+lwg)

    m = round(center_wavelength/2/dx)*dx
    # prob["margins"] = [[m, m], [m, m]]
    return prob


def port_number(port):
    s = str(port).split("@")[0]
    if s[0] == "o":
        s = s[1:]
    return int(s)


def mode_number(port):
    l = str(port).split("@")
    return 0 if len(l) == 1 else int(l[1])


def longname(s):
    o, i = s.split(",")
    po, pi = port_number(o), port_number(i)
    mo, mi = mode_number(o), mode_number(i)
    return f"o{po}@{mo},o{pi}@{mi}"
