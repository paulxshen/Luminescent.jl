import time
import pymeshfix
import json
from statistics import median
from .constants import *
import gdsfactory as gf
from gdsfactory.cross_section import Section
from .constants import *
from .materials import *
try:
    from IPython.display import display
except ImportError:
    pass
import pyvista as pv

# from .picToGDS import main

from math import ceil, cos, pi, sin, tan
import matplotlib.pyplot as plt
import numpy as np
from gdsfactory.generic_tech import LAYER_STACK, LAYER
import copy
import shutil
# from gdsfactory import LAYER_VIEWS

tol = .001


def arange(a, b, d):
    ret = np.linspace(a, b, round((b-a)/(d))+1).tolist()
    return ret


def trim(x, dx):
    return round(x/dx)*dx


def extend(endpoints, wm):
    v = endpoints[1]-endpoints[0]
    v = v/np.linalg.norm(v)
    return [
        (endpoints[0]-wm*v).tolist(), (endpoints[1]+wm*v).tolist()]


def portsides(c):
    ports = c.ports
    bbox = c.bbox_np()
    res = [[], [], [], []]
    xmin0, ymin0 = bbox[0]
    xmax0, ymax0 = bbox[1]
    for p in ports:
        x, y = np.array(p.center)/1e3

        if abs(x - xmin0) < tol:
            res[2].append(p.name)
        if abs(x - xmax0) < tol:
            res[0].append(p.name)
        if abs(y - ymin0) < tol:
            res[3].append(p.name)
        if abs(y - ymax0) < tol:
            res[1].append(p.name)
    return res


def add_bbox(c, layer, nonport_margin=0):  # , dx=None):
    bbox = c.bbox_np()
    xmin0, ymin0 = bbox[0]
    xmax0, ymax0 = bbox[1]
    l = xmax0-xmin0
    w = ymax0-ymin0

    # l = dx*np.ceil((xmax0-xmin0)/dx)
    # w = dx*np.ceil((ymax0-ymin0)/dx)

    # if dx is not None:
    #     if nonport_margin is None:
    #         nonport_margin = dx
    # if nonport_margin is None:
    #     nonport_margin = 0
    margin = nonport_margin
    xmin, ymin, xmax, ymax = xmin0-margin, ymin0 - \
        margin, xmin0+l+margin, ymin0+w+margin

    for p in c.ports:
        # p = c.ports[k]
        x, y = np.array(p.center)/1e3
        if abs(x - xmin0) < tol:
            xmin = x
        if abs(x - xmax0) < tol:
            xmax = x
        if abs(y - ymin0) < tol:
            ymin = y
        if abs(y - ymax0) < tol:
            ymax = y
    p = [(xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax)]
    _c = gf.Component()
    _c << c

    if type(layer[0]) is int:
        layer = [layer]
    for layer in layer:
        layer = tuple(layer)
        _c.add_polygon(p, layer=layer)
    for port in c.ports:
        _c.add_port(name=port.name, port=port)
    return _c

# def pic2gds(fileName, sizeOfTheCell, layerNum=1, isDither=False, scale=1):
    main(fileName, sizeOfTheCell, layerNum, isDither, scale)
    return "image.bmp", "image.gds"


def finish(c, name):
    c.add_label(name, position=c.bbox_np()[1])


def normal_from_orientation(orientation):
    return [cos(orientation/180*pi), sin(orientation/180*pi)]


def stl_to_array(mesh: pv.PolyData, dl: float, bbox):
    # x_min, x_max, y_min, y_max, z_min, z_max = mesh.bounds
    lb = [mesh.bounds[i] for i in [0, 2, 4]]
    ub = [mesh.bounds[i] for i in [1, 3, 5]]
    lims = [[a+dl/2+dl*ceil((p-a-dl/2)/dl), b-dl/2-dl*ceil((b-dl/2-q)/dl)]
            for (a, b, p, q) in zip(bbox[0], bbox[1], lb, ub)]
    xyz = [np.linspace(a, b, 1+round((b-a)/dl)) for (a, b) in lims]
    x, y, z = np.meshgrid(*xyz)

    # Create unstructured grid from the structured grid
    grid = pv.StructuredGrid(x, y, z)
    ugrid = pv.UnstructuredGrid(grid)

    # Get part of the mesh within the mesh's bounding surface.
    selection = ugrid.select_enclosed_points(
        mesh.extract_surface(),
        tolerance=0,
        check_surface=False,
    )
    mask = selection['SelectedPoints'].view(bool)
    mask = mask.reshape(x.shape, order='F')
    mask = np.array(mask)
    return mask


def material_voxelate(c, dl, zmin, zmax, layers, layer_stack, path):
    stacks = sum([[[v.mesh_order, v.material, tuple(layer), k]
                 for k, v in get_layers(layer_stack, layer, withkey=True)] for layer in layers], [])
    c.flatten()
    stacks = sorted(stacks, key=lambda x: -x[0])
    layer_stack_info = dict()
    # c.show()
    # raise NotImplementedError("This is a stub")
    lb, ub = c.bbox_np()
    # bbox = [[**lb, zmin], [**ub, zmax]]
    bbox = [[lb[0], lb[1], zmin], [ub[0], ub[1], zmax]]
    for i, stack in enumerate(stacks):
        m = stack[1]
        l1, l2 = layer = stack[2]
        k = stack[3]

        _layer_stack = copy.deepcopy(layer_stack)
        _layer_stack.layers.clear()

        d = copy.deepcopy(layer_stack.layers[k])
        if d.zmin <= zmax and d.bounds[1] >= zmin:
            _layer_stack.layers[k] = d
            d.zmin = max(zmin, d.zmin)
            d.thickness = min(zmax-d.zmin, d.thickness)
            # _d.bounds = (_d.zmin, _d.zmin+_d.thickness)
            # origin = (c.extract([layer]).bbox_np()-c.bbox_np())[0].tolist()
            STL = os.path.abspath(os.path.join(path, f"{k}.stl"))
            gf.export.to_stl(c, STL, layer_stack=_layer_stack)
            STL = os.path.join(path, f'{k}_{l1}_{l2}.stl')

            pymeshfix.clean_from_file(STL, STL)
            mesh = pv.read(STL)
            im = stl_to_array(mesh, dl, bbox)
            np.save(os.path.join(path, f'{k}.npy'), im)

            layer_stack_info[k] = {
                "layer": (l1, l2),
                "zmin": d.zmin,
                "thickness": d.thickness,
                "material": matname(m),
                "mesh_order": stack[0],
                # "origin": origin,
            }
    return layer_stack_info


def get_layers(layer_stack, layer, withkey=False):
    r = []
    for k, x in layer_stack.layers.items():
        l = x.layer
        if hasattr(l, "layer"):
            if tuple(l.layer) == tuple(layer):
                if withkey:
                    x = k, x
                # print(layer, k, x)
                r.append(x)
    if r:
        return r

    for k, x in layer_stack.layers.items():
        l = x.derived_layer
        if hasattr(l, "layer"):
            if tuple(l.layer) == tuple(layer):
                if withkey:
                    x = k, x
                # print(layer, k, x)
                r.append(x)
    return r


def wavelength_range(center, bandwidth, length=3):
    f1 = 1/(center+bandwidth/2)
    f2 = 1/(center-bandwidth/2)
    hw = (f2-f1)/2
    f1 = 1/center-hw
    f2 = 1/center+hw
    return sorted([1/x for x in np.linspace(f1, f2, length).tolist()])


def adjust_wavelengths(wavelengths, wl_res=.01):
    if len(wavelengths) == 1:
        return wavelengths, wavelengths[0], 1
    wavelengths = sorted(set(wavelengths))
    wavelengths0 = copy.deepcopy(wavelengths)
    center_wavelength = median(wavelengths)
    freqs0 = [center_wavelength/w for w in reversed(wavelengths)]
    nresfreq = round(1/min(np.diff([0]+freqs0)))

    while True:
        freqs = [round(f*nresfreq)/nresfreq for f in freqs0]
        wavelengths = [center_wavelength/f for f in reversed(freqs)]
        if max([abs(w-w0)/center_wavelength for w, w0 in zip(wavelengths, wavelengths0)]) <= wl_res:
            break
        nresfreq += 1
    print(
        f"wavelengths has been adjusted within `wl_res` to facilitate simulation:\n{wavelengths}")
    return wavelengths, center_wavelength, nresfreq


def save_prob(prob, path):
    path = os.path.abspath(path)
    bson_data = json.dumps(prob)
    # prob["component"] = c0

    path = prob["path"]
    if not os.path.exists(path):
        os.makedirs(path)
        #   compiling julia code...
        #   """)
    prob_path = os.path.join(path, "problem.json")
    print(prob_path)
    with open(prob_path, "w") as f:
        # Write the BSON data to the file
        f.write(bson_data)
    print("using simulation folder", path)
