# import dill
import math
from .constants import *
from .layers import *
from .utils import *
import bson
import gdsfactory as gf
from copy import deepcopy
# from time import time
import datetime
from math import cos, pi, sin
import os
import numpy as np

from sortedcontainers import SortedDict, SortedSet
from gdsfactory.generic_tech import LAYER_STACK, LAYER


def setup(c, study, dx, margin,
          bbox_layer=LAYER.WAFER,
          zmargin=None, zlims=None, core_layer=LAYER.WG,
          port_source_offset="auto", source_margin="auto",
          runs=[],  sources=[],
          layer_stack=LAYER_STACK, materials=MATERIALS,
          exclude_layers=[
              DESIGN_LAYER, GUESS], approx_2D=False, Courant=None,
          gpu=None, dtype=np.float32,
          plot=False, framerate=0,
          magic="", wd=os.path.join(os.getcwd(), "runs"), name=None, **kwargs):
    dx0 = dx
    # dx *= 2
    ratio = 4
    dl = dx/ratio
    dy = dx
    dz = dx
    prob = dict()
    prob = {**prob, **kwargs}
    prob["dx"] = dx
    prob["dy"] = dy
    prob["dz"] = dz
    prob["dl"] = dl
    prob["dx0"] = dx0
    prob["ratio"] = ratio
    prob["dtype"] = str(dtype)
    prob["timestamp"] = datetime.datetime.now().isoformat(
        timespec="seconds").replace(":", "-")
    prob["magic"] = magic
    prob["framerate"] = framerate
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

    # if zlims is None:
    # if core not in layer_stack.layers:
    #     raise ValueError("please set core=\"core_layer_name\"")
    # "\"core\" not found in layer_stack so `zlims` must be provided")
    d = get_layers(layer_stack, core_layer)[0]
    hcore = d.thickness
    zcore = d.zmin
    # hcore = round(d.thickness/dx)*dx
    # zlims = [zcore, zcore+hcore]
    thickness = hcore
    zcenter = zcore+hcore/2

    port_width = max([p.width/1e3 for p in c.ports])
    wg_margin = 1*port_width
    wg_zmargin = 1*thickness
    wmode = port_width+2*wg_margin
    hmode = thickness+2*wg_zmargin
    wmode = round(wmode / dx)*dx
    hmode = round(hmode / dz)*dz
    wg_margin = (wmode-port_width)/2
    wg_zmargin = (hmode-thickness)/2

    C = 1.4
    # if zmargin is None:
    zmargin = C*wg_zmargin
    # if margin is None:
    xmargin = ymargin = C * wg_margin

    h = thickness+2*zmargin
    h = round(h/dz)*dz
    zmargin = (h-thickness)/2
    zmin = zcore-zmargin

    prob["thickness"] = thickness
    prob["zcenter"] = zcenter
    prob["zmin"] = zmin
    prob["zcore"] = zcore

    # zlims = [zlims[0]-zmargin, zlims[1]+zmargin]
    # zcenter = zlims[0]+(zlims[1]-zlims[0])/2
    # zmin = zlims[0]

    length, width = c.bbox_np()[1]-c.bbox_np()[0]

    l = length+2*xmargin
    l = round(l/dx)*dx
    xmargin = (l-length)/2

    w = width+2*ymargin
    w = round(w/dx)*dx
    ymargin = (w-width)/2

    _c = gf.Component()
    _c << gf.components.bbox(
        component=c, left=xmargin, right=xmargin, top=ymargin, bottom=ymargin, layer=bbox_layer)
    for orientation in [0, 90, 180, 270]:
        if orientation in [0, 180]:
            length = xmargin
        else:
            length = ymargin
        c = gf.components.extend_ports(
            component=c, orientation=orientation, length=length)
    c0 = _c << c
    c = _c

    for run in runs:
        for port in run["sources"]:
            run["sources"][port]["center"] = (
                np.array(c0.ports[port].center)/1e3).tolist()
        # for port in run["monitors"]:
        #     x, y = c.ports[port].center
            # run["monitors"][port]["center"] = [x+xmargin, y+ymargin]

    prob["xmargin"] = xmargin
    prob["ymargin"] = ymargin
    prob["zmargin"] = zmargin

    # a = min([min([materials[d.material]["epsilon"] for d in get_layers(
    #     layer_stack, l)]) for l in bbox_layer])
    # b = max([materials[d.material]["epsilon"] for d in get_layers(
    #     layer_stack, core_layer)])
    # C = 4*math.sqrt(a/b)
    # if zmargin is None:
    #     zmargin = dx*round(C*thickness/dx)
    # port_width = max([p.width/1e3 for p in c.ports])
    # if margin is None:
    #     margin = trim(C*port_width, dx)
    layers = set(c.layers)-set(exclude_layers)

    # l, w = c.bbox_np()[1]-c.bbox_np()[0]
    center = [c.bbox_np()[0][0], c.bbox_np()[0][1]+w/2, zcenter]
    normal = [1, 0, 0]

    eps_3D = material_voxelate(c, dl, center, l, w, h,
                               normal, layers, layer_stack, materials)
    eps_2D = eps_3D[:, :, int(eps_3D.shape[2]/2)]
    prob["study"] = study

    if not name:
        l = [prob["timestamp"], study]
        name = "#".join(l)
    path = os.path.join(wd, name)
    prob["path"] = path
    prob["name"] = name

    prob["eps_3D"] = eps_3D.tolist()
    prob["eps_2D"] = eps_2D.tolist()
    prob["N"] = 2 if approx_2D else 3

    prob["hmode"] = hmode
    wmode = port_width+2*wg_margin
    # print(margin)
    neffmin = 1000000
    wavelengths = []
    # _c = add_bbox(c, layer=bbox_layer, nonport_margin=margin)
    for run in runs:
        for s in list(run["sources"].values())+list(run["monitors"].values()):
            s["mode_width"] = wmode
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
                        c, dl, center, wmode, hmode, normal, layers, layer_stack, materials)
                    # _modes, _modes1, neffs, _ = solve_modes(
                    _modes,  neffs = solve_modes(
                        eps, λ=wl, dx=dx, neigs=max(mode_numbers)+1, plot=plot)

                    for n in neffs:
                        if n < neffmin:
                            neffmin = n

                    modes = [{k: [np.real(mode[k].T).tolist(), np.imag(
                        mode[k].T).tolist()] for k in mode} for mode in _modes]
                    # modes1 = [{k: [np.real(mode[k]).tolist(), np.imag(
                    #     mode[k]).tolist()] for k in mode} for mode in _modes1]
                    # for mn in mode_numbers:
                    mode_solutions.append({
                        "modes": modes,
                        # "modes1": modes1,
                        "wavelength": wl,
                        # "mode_number": mn,
                        "ports": [s["port"]],
                        "size": [wmode, hmode],
                        "eps": eps.tolist(),
                        "width": s["width"],
                        "mode_width": wmode,
                        "zcenter": zcenter,
                    })
    wavelengths = sorted(set(wavelengths))
    wl = np.median(wavelengths)
    if port_source_offset == "auto":
        # port_source_offset = trim(2.6*wl/neffmin, dx)
        port_source_offset = trim(3*port_width, dx)
    prob["port_source_offset"] = port_source_offset
    if source_margin == "auto":
        source_margin = 2*dx
    prob["source_margin"] = source_margin

    bbox = c.bbox_np()
    source_ports = []
    nonsource_ports = []
    for p in c.ports:
        is_source = False
        for run in runs:
            for port in run["sources"]:
                if port == int(p.name[1]):
                    is_source = True
        if is_source:
            source_ports.append(p)
        else:
            nonsource_ports.append(p)

    prob["source_portsides"] = portsides(source_ports, bbox)
    prob["nonsource_portsides"] = portsides(nonsource_ports, bbox)

    prob["mode_solutions"] = mode_solutions
    prob["runs"] = runs
    prob["components"] = {
        "device": {
            "bbox": c.bbox_np().tolist(),
        }}

    prob["Courant"] = Courant
    if not os.path.exists(path):
        os.makedirs(path)
    # prob["wavelengths"] = wavelengths
    c.write_gds(os.path.join(path, "component.gds"))
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


def unpack_sparam_key(k):
    o, i = k.split(",")
    po, pi = port_number(o), port_number(i)
    mo, mi = mode_number(o), mode_number(i)
    return f"o{po}", mo, f"o{pi}", mi


def long_sparam_key(k):
    po, mo, pi, mi = unpack_sparam_key(k)
    return f"o{po}@{mo},o{pi}@{mi}"
