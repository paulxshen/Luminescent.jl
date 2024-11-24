from .constants import *
from .layers import *
from .utils import *
import json
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
          zmargin2=None, zlims=None, core_layer=LAYER.WG,
          port_source_offset="auto", source_margin="auto",
          runs=[],  sources=[],
          layer_stack=LAYER_STACK, materials=MATERIALS,
          exclude_layers=[
              DESIGN_LAYER, GUESS], N=3, Courant=None,
          gpu=None, dtype=np.float32,
          plot=False, framerate=0,
          magic="", wd=os.path.join(os.getcwd(), "runs"), name=None, **kwargs):
    prob = dict()
    dx0 = dx
    # dx *= 2
    ratio = 6
    dl = dx/ratio
    # dl = .01
    # ratio = int(dx/dl)
    dy = dx
    dz = 1 * dx

    if not name:
        l = [prob["timestamp"], study]
        name = "#".join(l)
    path = os.path.join(wd, name)
    prob["path"] = path
    prob["name"] = name

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

    d = get_layers(layer_stack, core_layer)[0]
    hcore = d.thickness
    zcore = d.zmin

    zmargin1 = 2*hcore
    hmode = hcore+2*zmargin1
    hmode = trim(hmode, 2*dz)
    zmargin1 = (hmode-hcore)/2

    zmode = zcore-zmargin1
    zcenter = zcore+hcore/2

    zmargin2 = 2.5*hcore
    zmargin2 = trim(zmargin2,  dz)
    h = hcore+2*(zmargin1+zmargin2)
    zmin = zcore-zmargin2-zmargin1
    zmax = zmin+h

    source_ports = []
    nonsource_ports = []
    for p in c.ports:
        is_source = False
        for run in runs:
            for port in run["sources"]:
                if port == p.name:
                    is_source = True
        if is_source:
            source_ports.append(p.name)
        else:
            nonsource_ports.append(p.name)
#
    port_width = max([p.width/1e3 for p in c.ports])
    ps = portsides(c)
    xmargin = ymargin = 2*port_width
    source_port_margin = 6*port_width
    margins = []
    for p in ps:
        if set(p).intersection(source_ports):
            margins.append(source_port_margin)
        # elif p:
            # margins.append(border_margin)
        else:
            margins.append(xmargin)
    l0, w0 = c.bbox_np()[1]-c.bbox_np()[0]

    _l = l0+margins[0]+margins[2]
    l = trim(_l, 2*dx)
    margins[0] += (l-_l)

    _w = w0+margins[1]+margins[3]
    w = trim(_w, 2*dy)
    margins[1] += (w-_w)

    #
    modemargin = 1*port_width
    wmode = port_width+2*modemargin
    wmode = trim(wmode, 2*dx)
    modemargin = (wmode-port_width)/2

    prob["hcore"] = hcore
    prob["zcenter"] = zcenter
    prob["zmin"] = zmin
    prob["zcore"] = zcore
    prob["xmargin"] = xmargin
    prob["ymargin"] = ymargin
    prob["zmargin2"] = zmargin2
    prob["L"] = [l, w, h]

    _c = gf.Component()
    kwargs = dict()
    for (orientation, side, length, p) in zip([0, 90, 180, 270], ["right", "top", "left", "bottom"], margins, ps):
        if p:
            c = gf.components.extend_ports(
                c, None, length, orientation=orientation)
        else:
            kwargs[side] = length
    _c << gf.components.bbox(component=c, layer=bbox_layer, **kwargs)
    c0 = _c << c
    c = _c
    # c.plot()
    # c.show()

    lb = c.bbox_np()[0].tolist()
    center = [c.bbox_np()[0][0], c.bbox_np()[0][1]+w/2, zcenter]

    xs = arange(lb[0], lb[0]+l, dx)
    ys = arange(lb[1], lb[1]+w, dy)
    # zs = sorted(list(set(
    #     arange(zmin, zcore-zmargin1, 4*dz) +
    #     arange(zcore-zmargin1, zcore-2*dz, 2*dz) +
    #     arange(zcore-2*dz, zcore+hcore+2*dz, dz) +
    #     arange(zcore+hcore+2*dz, zcore+hcore+zmargin1, 2*dz) +
    #     arange(zcore+hcore+zmargin1, zmax, 4*dz)
    # )))
    z0 = zmin
    z1 = zmode
    z2 = zmode+hmode
    z3 = zmax
    zs = sorted(list(set(
        arange(z0, z1, 2*dz) +
        arange(z1, z2, dz) +
        arange(z2, z3, 2*dz)
    )))
    prob["xs"] = xs
    prob["ys"] = ys
    prob["zs"] = zs

    for run in runs:
        for port in run["sources"]:
            run["sources"][port]["center"] = (
                np.array(c0.ports[port].center)/1e3).tolist()

    # a = min([min([materials[d.material]["epsilon"] for d in get_layers(
    #     layer_stack, l)]) for l in bbox_layer])
    # b = max([materials[d.material]["epsilon"] for d in get_layers(
    #     layer_stack, core_layer)])
    # C = 4*math.sqrt(a/b)
    # if zmargin2 is None:
    #     zmargin2 = dx*round(C*hcore/dx)
    # port_width = max([p.width/1e3 for p in c.ports])
    # if margin is None:
    #     margin = trim(C*port_width, dx)
    layers = set(c.layers)-set(exclude_layers)

    normal = [1, 0, 0]

    temp = os.path.join(path, "temp")
    os.makedirs(temp, exist_ok=True)
    layer_stack_info = material_voxelate(
        c, dl, zmin, zmax, layers, layer_stack, temp)
    prob["layer_stack"] = layer_stack_info
    prob["materials"] = materials
    prob["study"] = study

    prob["N"] = N

    prob["hmode"] = hmode
    prob["wmode"] = wmode
    prob["zmode"] = zmode
    wmode = port_width+2*modemargin
    # print(margin)
    neffmin = 1000000
    wavelengths = []
    # _c = add_bbox(c, layer=bbox_layer, nonport_margin=margin)
    for run in runs:
        for s in list(run["sources"].values())+list(run["monitors"].values()):
            s["mode_width"] = wmode
            for wl in s["wavelength_mode_numbers"]:
                wavelengths.append(wl)
    wavelengths = sorted(set(wavelengths))
    wl = np.median(wavelengths)
    if source_margin == "auto":
        source_margin = 2*dx
    prob["source_margin"] = source_margin

    bbox = c.bbox_np()

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


def port_name(port):
    s = str(port).split("@")[0]
    if s[0] == "o":
        return s
    return f"o{s}"


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
    po, pi = port_name(o), port_name(i)
    mo, mi = mode_number(o), mode_number(i)
    return po, mo, pi, mi


def long_sparam_key(k):
    po, mo, pi, mi = unpack_sparam_key(k)
    return f"{po}@{mo},{pi}@{mi}"
