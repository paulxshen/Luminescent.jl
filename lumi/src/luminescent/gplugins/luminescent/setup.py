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


def setup(c, study,   dx,
          margin,  zmargin, port_source_offset="auto", source_margin="auto", name="",
          runs=[],  sources=[], layer_stack=LAYER_STACK, layer_views=LAYER_VIEWS, exclude_layers=[
              LAYER.DESIGN, GUESS], approx_2D=False, Courant=None,
          gpu=None, dtype=np.float32,
          path=PATH, plot=False, **kwargs):
    prob = dict()
    prob = {**prob, **kwargs}
    prob["dtype"] = str(dtype)
    prob["timestamp"] = datetime.datetime.now().isoformat(
        timespec="seconds").replace(":", "-")
    prob["name"] = name
    prob["gpu_backend"] = gpu if gpu else ""
    ports = {
        p.name: {
            "center": p.center,
            "normal": [cos(p.orientation/180*pi), sin(p.orientation/180*pi)],
        }
        for p in c.get_ports_list(prefix="o")
    }
    prob["ports"] = ports

    mode_solutions = []
    layers = set(c.layers)-set(exclude_layers)
    layers = sorted(layers, key=lambda layer: -get_layer(
        layer_stack, layer).mesh_order)

    hcore = layer_stack.layers["core"].thickness
    # hcore = round(layer_stack.layers["core"].thickness/dx)*dx
    if zmargin is None:
        zmargin = dx*round(1.5*hcore/dx)
    zcenter = layer_stack.layers["core"].zmin+hcore/2
    h = hcore+2*zmargin
    prob["mode_height"] = h
    l, w = c.bbox_np()[1]-c.bbox_np()[0]
    # l, w = round(l/dx)*dx, round(w/dx)*dx
    center = [c.bbox_np()[0][0], c.bbox_np()[0][1]+w/2, zcenter]
    normal = [1, 0, 0]
    eps = material_voxelate(c, dx, center, l, w, h,
                            normal, layers, layer_stack)
    eps_2D = eps[:, :, int(eps.shape[2]/2)]
    prob["study"] = study
    prob["path"] = os.path.join(
        path, prob["study"] + "_" + prob["timestamp"])

    prob["eps_3D"] = eps.tolist()
    prob["eps_2D"] = eps_2D.tolist()
    prob["zmin"] = -zmargin
    prob["d"] = 2 if approx_2D else 3

    neffmin = 1000000
    wavelengths = []
    for run in runs:
        for s in list(run["sources"].values())+list(run["monitors"].values()):
            for wl in s["wavelength_mode_numbers"]:
                wavelengths.append(wl)
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
                    _modes, _modes1, neffs, _ = solve_modes(
                        eps, λ=wl, dx=dx, neigs=max(mode_numbers)+1, plot=plot)

                    for n in neffs:
                        if n < neffmin:
                            neffmin = n

                    modes = [{k: [np.real(mode[k].T).tolist(), np.imag(
                        mode[k].T).tolist()] for k in mode} for mode in _modes]
                    modes1 = [{k: [np.real(mode[k]).tolist(), np.imag(
                        mode[k]).tolist()] for k in mode} for mode in _modes1]
                    # for mn in mode_numbers:
                    mode_solutions.append({
                        "modes": modes,
                        "modes1": modes1,
                        "wavelength": wl,
                        # "mode_number": mn,
                        "ports": [s["port"]],
                        "size": [w, h],
                        "eps": eps.tolist(),
                        "width": w,
                        "zcenter": zcenter,
                    })
    wavelengths = sorted(set(wavelengths))
    wl = np.median(wavelengths)
    if port_source_offset == "auto":
        port_source_offset = 2*wl/neffmin
    prob["port_source_offset"] = port_source_offset
    if source_margin == "auto":
        source_margin = .2*wl/neffmin
    prob["source_margin"] = source_margin

    prob["margin"] = margin
    prob["zmargin"] = zmargin
    prob["portsides"] = portsides(c)

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
    prob["Courant"] = Courant
    prob["component"] = c
    # λc = (wavelengths[0]+wavelengths[-1])/2
    # prob["wavelengths"] = wavelengths
    prob["mode_solutions"] = mode_solutions
    # prob["path_length"] = 2*(l+r+lwg)
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
