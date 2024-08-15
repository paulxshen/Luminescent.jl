import dill
import math
from .constants import *
from .layers import *
from .utils import *
import bson
import gdsfactory as gf
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
from gdsfactory.generic_tech import LAYER_STACK, LAYER


def setup(c, study,   dx, margin,
          bbox_layer=LAYER.WAFER,
          zmargin=None, zlims=None, core_layer=LAYER.WG,
          port_source_offset="auto", source_margin="auto", name="",
          runs=[],  sources=[],
          layer_stack=LAYER_STACK, materials=MATERIALS,
          exclude_layers=[
              DESIGN_LAYER, GUESS], approx_2D=False, Courant=None,
          gpu=None, dtype=np.float32,
          path=RUNS_PATH, plot=False, **kwargs):
    if type(bbox_layer[0]) is int:
        bbox_layer = (bbox_layer,)
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

    if zlims is None:
        # if core not in layer_stack.layers:
        #     raise ValueError("please set core=\"core_layer_name\"")
        # "\"core\" not found in layer_stack so `zlims` must be provided")
        d = get_layers(layer_stack, core_layer)[0]
        hcore = d.thickness
        zmin = d.zmin
        # hcore = round(d.thickness/dx)*dx
        zlims = [zmin, zmin+hcore]
    thickness = zlims[1]-zlims[0]

    a = min([min([materials[d.material].epsilon for d in get_layers(
        layer_stack, l)]) for l in bbox_layer])
    b = max([materials[d.material].epsilon for d in get_layers(
        layer_stack, core_layer)])
    C = 5*math.sqrt(a/b)
    if zmargin is None:
        zmargin = dx*round(C*thickness/dx)
    port_width = max([p.width/1e3 for p in c.ports])
    if margin is None:
        margin = trim(C*port_width, dx)

    c = add_bbox(c, layer=bbox_layer, nonport_margin=margin)
    layers = set(c.layers)-set(exclude_layers)

    # zlims = [zlims[0]-dx, zlims[0]+np.ceil(thickness/dx)*dx+dx]
    zlims = [zlims[0]-zmargin, zlims[1]+zmargin]
    h = zlims[1]-zlims[0]
    zcenter = zlims[0]+(zlims[1]-zlims[0])/2
    zmin = zlims[0]

    l, w = c.bbox_np()[1]-c.bbox_np()[0]
    # l, w = round(l/dx)*dx, round(w/dx)*dx
    center = [c.bbox_np()[0][0], c.bbox_np()[0][1]+w/2, zcenter]
    normal = [1, 0, 0]

    eps = material_voxelate(c, dx, center, l, w, h,
                            normal, layers, layer_stack, materials)
    eps_2D = eps[:, :, int(eps.shape[2]/2)]
    prob["study"] = study

    l = [prob["timestamp"], study]
    if name:
        l.append(name)
    path = os.path.join(path, "#".join(l))
    prob["path"] = path

    prob["eps_3D"] = eps.tolist()
    prob["eps"] = prob["eps_3D"]
    prob["eps_2D"] = eps_2D.tolist()
    prob["zmin"] = zmin
    prob["d"] = 2 if approx_2D else 3

    prob["mode_height"] = h
    w = port_width+2*margin
    neffmin = 1000000
    wavelengths = []
    # _c = add_bbox(c, layer=bbox_layer, nonport_margin=margin)
    for run in runs:
        for s in list(run["sources"].values())+list(run["monitors"].values()):
            s["mode_width"] = w
            for wl in s["wavelength_mode_numbers"]:
                wavelengths.append(wl)
                duplicate = False
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
                        c, dx, center, w, h, normal, layers, layer_stack, materials)
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
                        "width": s["width"],
                        "mode_width": w,
                        "zcenter": zcenter,
                    })
    wavelengths = sorted(set(wavelengths))
    wl = np.median(wavelengths)
    if port_source_offset == "auto":
        port_source_offset = trim(2*wl/neffmin, dx)
    prob["port_source_offset"] = port_source_offset
    if source_margin == "auto":
        source_margin = 2*dx
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
    # prob["component"] = c
    dill.dump(c, os.path.join(path, "comp.pk"))
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
