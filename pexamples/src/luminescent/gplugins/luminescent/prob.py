import subprocess
import sys
from altair import layer
from networkx import center
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
from sortedcontainers import SortedDict, SortedSet
from .generic_tech import LAYER_STACK, LAYER_MAP, LAYER_VIEWS
import gdsfactory as gf
import bson
import scipy as sp
from .utils import *
from .layers import *


# include_layers=[layer_map.WG, layer_map.WAFER], ** kwargs):


def setup(c, study, center_wavelength,  dx,
          waveguide_height_margin,
          waveguide_width_margin, name="",
          runs=[], wavelengths=[], sources=[], layer_stack=LAYER_STACK, layer_views=LAYER_VIEWS, exclude_layers=[LAYER_MAP.DESIGN, GUESS], approx_2D=False):
    prob = dict()
    prob["name"] = name
    ports = {
        k: {
            "center": p.center.flatten().tolist(),
            #  "lb": p.endpoints[0].flatten().tolist(),
            #  "ub": p.endpoints[1].flatten().tolist(),
            "normal": [cos(p.orientation/180*pi), sin(p.orientation/180*pi)],
            "endpoints": p.endpoints.tolist(),
            # **monitors[k],
        }
        for k, p in zip(c.ports, c.ports.values(), ) if k.startswith("o")
    }
    prob["ports"] = ports

    mode_solutions = []
    if study == "sparams":
        1
    else:
        1
        # prob["sources"] = c.metadata["sources"]
        # prob["monitors"] = c.metadata["monitors"]
    layers = set(c.layers)-set(exclude_layers)
    layers = sorted(layers, key=lambda layer: next(
        x for x in layer_stack.layers.values() if x.layer == layer).mesh_order)

    hcore = layer_stack.layers["core"].thickness
    zcenter = layer_stack.layers["core"].zmin+hcore/2
    h = hcore+2*waveguide_height_margin
    prob["mode_height"] = h
    l, w = c.bbox[1]-c.bbox[0]
    center = [c.bbox[0][0], c.bbox[0][1]+w/2, zcenter]
    normal = [1, 0, 0]
    eps = material_voxelate(c, dx, center, l, w, h,
                            normal, layers, layer_stack)
    eps_2D = eps[:, :, int(eps.shape[2]/2)]
    prob["eps_3D"] = eps.tolist()
    prob["eps_2D"] = eps_2D.tolist()

    for run in runs:
        for s in list(run["sources"].values())+list(run["monitors"].values()):
            for wl in s["wavelength_mode_numbers"]:
                for mn in s["wavelength_mode_numbers"][wl]:
                    duplicate = False
                    wwg = s["width"]
                    for m in mode_solutions:
                        if m["wavelength"] == wl and m["mode_number"] == mn and m["width"] == s["width"]:
                            m["ports"].append(s["port"])
                            duplicate = True
                    if not duplicate:

                        w = wwg+2*waveguide_width_margin
                        normal = s["normal"]+[0]
                        center = s["center"]+[zcenter]

                        center = np.array(center)-.001*np.array(normal)
                        eps = material_slice(
                            c, dx, center, w, h, normal, layers, layer_stack, layer_views=layer_views)
                        _modes = solve_modes(eps, λ=wl, dx=dx,)
                        # mode = _modes[0]
                        modes = [{k: [np.real(mode[k]).tolist(), np.imag(
                            mode[k]).tolist()] for k in mode} for mode in _modes]
                        mode_solutions.append({
                            "modes": modes,
                            "wavelength": wl,
                            "mode_number": mn,
                            "ports": [s["port"]],
                            "size": [w, h],
                            "eps": eps.tolist(),
                            "width": wwg,
                        })
    prob["mode_solutions"] = mode_solutions
    prob["runs"] = runs

    device_svg = write_img("device", c,
                           hidden_layer=(DESIGN, GUESS))
    prob["components"] = {
        "device": {
            "svg": device_svg,
            "bbox": c.bbox.tolist(),
        }}
    prob["dx"] = dx
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


def sparams_problem(c, waveguide_width_margin=.4, waveguide_height_margin=.2, dx=.05, wavelengths=[1.55], center_wavelength=None, keys=None,
                    approx_2D=False, layer_stack=LAYER_STACK, **kwargs):
    d = 2 if approx_2D else 3
    prob = dict()

    ports = [int(p[1])
             for p in sorted(filter(lambda x: x[0] == "o", c.ports.keys()))]
    if keys is None:
        io = {i: ports for i in ports}
    else:
        io = SortedDict()
        for k in keys:
            (o, i) = longname(k).split(",")
            if i not in io:
                io[i] = SortedSet()
            io[i].add(o)
        for i in io:
            if i not in io[i]:
                io[i].add(i)
    runs = [{
        "d": d,
        "sources": {
            port_number(i): {
                "center": c.ports[f"o{port_number(i)}"].center.tolist(),
                "width": c.ports[f"o{port_number(i)}"].width,
                "normal": normal_from_orientation(c.ports[f"o{port_number(i)}"].orientation),
                "endpoints": extend(c.ports[f"o{port_number(i)}"].endpoints, waveguide_width_margin),
                "wavelength_mode_numbers": {w: [mode_number(i)] for w in wavelengths},
                "port": port_number(i),
            }},
        "monitors": {
            port_number(o): {
                "port": port_number(o),
                "normal": normal_from_orientation(c.ports[f"o{port_number(o)}"].orientation),
                "center": c.ports[f"o{port_number(o)}"].center.tolist(),
                "width": c.ports[f"o{port_number(o)}"].width,
                "endpoints": extend(c.ports[f"o{port_number(o)}"].endpoints, waveguide_width_margin),
                "wavelength_mode_numbers": {w: [mode_number(o)] for w in wavelengths},
            } for o in io[i]},
    } for i in io]

    prob = setup(c, study="sparams", center_wavelength=center_wavelength, dx=dx,
                 runs=runs, waveguide_width_margin=waveguide_width_margin,
                 waveguide_height_margin=waveguide_height_margin, layer_stack=layer_stack, approx_2D=approx_2D, **kwargs)
    # prob["runs"] = runs
    prob["study"] = "sparams"
    return prob


def inverse_design_problem(c, sparam_targets, lmin=.1, symmetries=[],
                           maxiters=25,
                           design_region_layer=LAYER_MAP.DESIGN,
                           design_guess_layer=LAYER_MAP.GUESS,
                           design_layer=LAYER_MAP.WG,
                           void_layer=LAYER_MAP.WGCLAD,
                           layer_stack=LAYER_STACK, **kwargs):
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

    prob["designs"] = [
        {
            "layer": design_layer,
            "bbox": [np.min(p, 0).tolist(), np.max(p, 0).tolist()],
            "symmetries": s,
            "guess": None,
            "lmin": lmin,
        } for p, s in zip(polys, symmetries)
    ]

    prob["design_config"] = dict()
    # prob[""]
    l = next(x for x in layer_stack.layers.values() if x.layer == design_layer)
    d = {"thickness": l.thickness, "material": l.material, "zmin": l.zmin}
    # d = vars(prob["design_layer"])
    d["ϵ"] = MATERIAL_EPS[d["material"]]
    prob["design_config"]["fill"] = d

    l = next(x for x in layer_stack.layers.values() if x.layer == void_layer)
    d = {"thickness": l.thickness, "material": l.material, "zmin": l.zmin}
    d["ϵ"] = MATERIAL_EPS[d["material"]]
    prob["design_config"]["void"] = d

    prob["design_layer"] = d
    prob["targets"] = targets
    prob["maxiters"] = maxiters
    return prob


def write_sparams(*args, **kwargs):
    return solve(sparams_problem(*args, **kwargs))


def solve(prob,):
    bson_data = bson.dumps(prob)
    dir = "_temp"
    if not os.path.exists(dir):
        os.makedirs(dir)
    prob_path = os.path.join(dir, "_prob.bson")
    with open(prob_path, "wb") as f:
        # Write the BSON data to the file
        f.write(bson_data)
    subprocess.run(
        [f"julia", os.path.join(os.path.dirname(os.path.abspath(__file__)), "run.jl"), prob_path,])
    sol = bson.loads(open(os.path.join(dir, "_sol.bson"), "rb").read())["sol"]
    if prob["study"] == "sparams":
        sol["sparams"] = {k: (v[0]+1j*v[1])
                          for k, v in sol["sparams"].items()}
    elif prob["study"] == "inverse_design":
        sol["designs"] = [np.array(d) for d in sol["designs"]]
    return sol
    # @save "$(@__DIR__)/layout.bson" device_mask sources ports designs dx λ ϵbase ϵclad ϵcore hbase hwg hclad modes l w
# if __module__
# main(1, .5, .22)


def apply_design(c, designs):
    res = gf.Component()
    c = res << c
    for d in designs:
        layer = d["layer"]
        x0, y0 = d["bbox"][0]
        x1, y1 = d["bbox"][1]
        p = res.add_polygon(
            [(x0, y0), (x1, y0), (x1, y1), (x0, y1)], layer=layer)
        res = gf.boolean(c, p, "not")
    return res
