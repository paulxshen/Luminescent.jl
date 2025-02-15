import platform
import subprocess
from ..constants import *
from ..layers import *
from ..utils import *
from ..materials import *
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


def setup(path, c, study, nres, center_wavelength,
          bbox_layer=BBOX_LAYER,
          zmargin2=None, zlims=None, core_layer=LAYER.WG,
          port_source_offset="auto", port_margin="auto",
          runs=[],  sources=[],
          layer_stack=SOI, materials=dict(),
          default_material='SiO2',
          exclude_layers=[
              DESIGN_LAYER, GUESS], Courant=None,
          gpu=None, dtype=np.float32,
          plot=False, framerate=None,
          magic="", wd=os.path.join(os.getcwd(), "runs"), name=None,
          Ttrans=None,
          approx_2D_mode=False):
    RATIO = 4

    materials = {**MATERIALS, **materials}
    prob = dict()
    if approx_2D_mode:
        N = 2
        prob["approx_2D_mode"] = approx_2D_mode
    else:
        N = 3
        prob["approx_2D_mode"] = None
    dy = dx = center_wavelength/nres
    dl = dx/RATIO
    dz = 1 * dx

    prob["class"] = "pic"
    prob["Ttrans"] = Ttrans
    prob["path"] = path
    prob["name"] = name
    prob["center_wavelength"] = center_wavelength
    prob["dx"] = dx
    prob["dy"] = dy
    prob["dz"] = dz
    prob["dl"] = dl
    # ratio =
    # prob["ratio"] = ratio
    prob["dtype"] = str(dtype)
    prob["timestamp"] = datetime.datetime.now().isoformat(
        timespec="seconds").replace(":", "-")
    prob["magic"] = magic
    prob["framerate"] = framerate

    gpu_backend = gpu
    # if gpu_backend:s
    prob["gpu_backend"] = gpu_backend

    ports = {
        p.name: {
            "center": (np.array(p.center)/1e3).tolist(),
            "normal": [cos(p.orientation/180*pi), sin(p.orientation/180*pi)],
            "tangent": [-sin(p.orientation/180*pi), cos(p.orientation/180*pi)],
            'width': p.width/1e3,
        }
        for p in c.get_ports_list(prefix="o")
    }
    for (k, v) in ports.items():
        z, n = [0, 0, 1], [*v['normal'], 0]
        t = np.cross(z, n).tolist()
        v['frame'] = [t, z, n]
        v["dimensions"] = np.abs(v['width']*np.array(v['tangent'])).tolist(),
    prob["ports"] = ports
    mode_solutions = []

    d = get_layers(layer_stack, core_layer)[0]
    hcore = d.thickness
    zcore = d.zmin

    zmargin1 = 3*hcore
    hmode = hcore+2*zmargin1
    hmode = trim(hmode, 2*dz)
    zmargin1 = (hmode-hcore)/2

    zmode = zcore-zmargin1
    zcenter = zcore+hcore/2

    zmargin2 = 2*hcore
    zmargin2 = trim(zmargin2,  2*dz)
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

    source_port_margin = 3 * port_width if N == 2 else 6*port_width

    port_margin = 2*dx
    margins = []
    for p in ps:
        if set(p).intersection(source_ports):
            margins.append(source_port_margin+port_margin)
        elif set(p).intersection(nonsource_ports):
            margins.append(port_margin)
        else:
            assert not p
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
    prob["zmax"] = zmax
    prob["zcore"] = zcore
    prob["xmargin"] = xmargin
    prob["ymargin"] = ymargin
    prob["zmargin2"] = zmargin2
    # prob["L"] = [l, w, h]

    _c = gf.Component()
    kwargs = dict()
    for (orientation, side, length, p) in zip([0, 90, 180, 270], ["right", "top", "left", "bottom"], margins, ps):
        if p:
            c = gf.components.extend_ports(
                c, None, length, orientation=orientation)
        else:
            kwargs[side] = length
    _c << gf.components.bbox(component=c, layer=bbox_layer, **kwargs)
    # _c = gf.Component()
    _c << c
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

    layers = set(c.layers)-set(exclude_layers)

    TEMP = os.path.join(path, "TEMP")
    os.makedirs(TEMP, exist_ok=True)
    SURFACES = os.path.join(path, 'surfaces')
    os.makedirs(SURFACES, exist_ok=True)

    layer_stack_info = material_voxelate(
        c, dl, zmin, zmax, layers, layer_stack, SURFACES)
    dir = os.path.dirname(os.path.realpath(__file__))
    print(dir)
    fn = os.path.join(dir, "solvemodes.py")
    assert os.path.exists(fn)
    assert os.path.exists(TEMP)
    if platform.system() == "Windows":
        os.system(f"copy /Y {fn} {TEMP}")
    else:
        subprocess.run(["cp", fn, TEMP])
    prob["layer_stack"] = layer_stack_info
    # materials = set([v["material"] for v in layer_stack_info.values()])
    # d = {
    #     k: materials[k] for k in materials
    # }
    prob["study"] = study
    # mateps = {k: np.real(v) for k, v in d.items()}
    # materials = {
    #     k: {} for k in materials
    # }
    # for k, v in mateps.items():
    #     materials[k]["epsilon"] = v
    prob["materials"] = materials

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
        for k, v in list(run["sources"].items())+list(run["monitors"].items()):
            p = ports[k]
            v['dimensions'] = p['dimensions']
            v['frame'] = p['frame']
            for center_wavelength in v["wavelength_mode_numbers"]:
                wavelengths.append(center_wavelength)

            if k in run["sources"]:
                ct = np.array(p['center'])
                n = np.array(p['normal'])
                v["center"] = (ct-n*port_margin).tolist()
            else:
                v['center'] = p['center']

    wavelengths = sorted(set(wavelengths))

    prob["mode_solutions"] = mode_solutions
    prob["runs"] = runs
    prob["components"] = {
        "device": {
            "bbox": c.bbox_np().tolist(),
        }}
    bbox = c.bbox_np().tolist()

    prob['bbox'] = bbox
    # prob['epdefault'] = materials[layer_stack['default']['material']]['epsilon']
    prob['epdefault'] = materials[default_material]['epsilon']

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
