import stltovoxel
import copy
from .constants import *
import gdsfactory as gf
from gdsfactory.cross_section import Section
from .constants import *
from .materials import MATERIALS
try:
    from IPython.display import display
except ImportError:
    pass

from re import T
from PIL import Image
import imageio.v3 as iio
import textwrap
from cairosvg import svg2png
from regex import F
from .picToGDS import main
import rasterio.features
import shapely
from shapely import GeometryCollection, Polygon
from shapely.ops import clip_by_rect
from shapely.plotting import plot_polygon, plot_points
from .picToGDS import main
import io
import os
import cv2
import gdstk
import pylab
import EMpy
# import EMpy as em
import numpy
from math import cos, pi, sin, tan
import trimesh
import matplotlib.pyplot as plt
import numpy as np
from gdsfactory.generic_tech import LAYER_STACK, LAYER
# from gdsfactory import LAYER_VIEWS

tol = .001


def trim(x, dx):
    return round(x/dx)*dx


def extend(endpoints, wm):
    v = endpoints[1]-endpoints[0]
    v = v/np.linalg.norm(v)
    return [
        (endpoints[0]-wm*v).tolist(), (endpoints[1]+wm*v).tolist()]


def portsides(ports, bbox):
    res = [[False, False], [False, False]]
    xmin0, ymin0 = bbox[0]
    xmax0, ymax0 = bbox[1]
    for p in ports:
        x, y = np.array(p.center)/1e3
        if abs(x - xmin0) < tol:
            res[0][0] = True
        if abs(x - xmax0) < tol:
            res[1][0] = True
        if abs(y - ymin0) < tol:
            res[0][1] = True
        if abs(y - ymax0) < tol:
            res[1][1] = True
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


def solve_modes(eps, λ, dx, neigs=1, plot=False):
    tol = 1e-4
    m, n = eps.shape
    # print(m, n)
    x = numpy.linspace(0.5*dx, (m-.5)*dx, m)
    y = numpy.linspace(0.5*dx, (n-.5)*dx, n)

    i = round(n/2)
    y1 = y[i:i+4]
    e = eps[:, i]
    eps1 = np.stack([e, e, e, e]).T

    def ϵfunc(x_, y_):
        m, n = len(x_), len(y_)
        return cv2.resize(eps, dsize=(n, m), interpolation=cv2.INTER_CUBIC)

    def ϵfunc1(x_, y_):
        m, n = len(x_), len(y_)
        return cv2.resize(eps1, dsize=(n, m), interpolation=cv2.INTER_CUBIC)
        # return ϵfunc(x_, y_)[i:i+len(y_), :]

    tol = 1e-6
    # plot = True
    # neigs = 2
    solver = EMpy.modesolvers.FD.VFDModeSolver(
        λ, x, y, ϵfunc,  "0000").solve(neigs, tol)
    solver1D = EMpy.modesolvers.FD.VFDModeSolver(
        λ, x, y1, ϵfunc1,  "AA00").solve(2*neigs, tol)
    # solvers = EMpy.modesolvers.FD.VFDModeSolver(
    #     λ, x, y1, ϵfunc1, "SS00").solve(neigs, tol)

    # fig = pylab.figure()
    # fig.add_subplot(2, 3, 1)
    # Hy = numpy.transpose(solver.modes[-1].get_field("Hy", x, y))
    # pylab.imshow(abs(Hy))
    # fig.add_subplot(2, 3, 2)
    # Hy = numpy.transpose(solver1D.modes[-1].get_field("Hy", x, y1))
    # pylab.imshow(abs(Hy))
    # fig.add_subplot(2, 3, 4)
    # e = numpy.transpose(ϵfunc(x, y))
    # pylab.imshow(e)
    # fig.add_subplot(2, 3, 5)
    # e = numpy.transpose(ϵfunc1(x, y1))
    # pylab.imshow(e)
    # pylab.show()

    modes = [{k: m.get_field(k, x, y) for k in [
        "Ex", "Ey", "Ez", "Hx", "Hy", "Hz"]} for m in solver.modes]
    neffs = [np.real(m.neff) for m in solver.modes]

    def f(m):
        a = m.get_field("Hy", x, y1)
        a = abs(a)
        return (max(np.std(a, 1)) / np.max(a)) < .1
    v = sorted(solver1D.modes, key=lambda x: -np.abs(x.neff))
    v = filter(f, v)
    modes1D = [{k: m.get_field(k, x, y1)[:, 0] for k in [
        "Ex", "Ey", "Ez", "Hx", "Hy", "Hz"]} for m in v]
    neffs1D = [np.real(m.neff) for m in v]
    # "Ex", "Ey", "Ez", "Hx", "Hy", "Hz"]} for m in sorted(solver1D.modes+solvers.modes, key=lambda x: -np.abs(x.neff))]
    # print(solver.modes[0].get_field("Ex", x, y).shape)

    return modes, neffs, modes1D, neffs1D


def s2svg(area, bbox, dx):
    # specify margin in coordinate units

    width = bbox[2] - bbox[0]
    height = bbox[3] - bbox[1]
    props = {
        'version': '1.1',
        'baseProfile': 'full',
        'width': '{width:.0f}px'.format(width=round(width/dx)),
        'height': '{height:.0f}px'.format(height=round(height/dx)),
        'viewBox': '%.1f,%.1f,%.1f,%.1f' % (bbox[0], bbox[1], width, height),
        'xmlns': 'http://www.w3.org/2000/svg',
        'xmlns:ev': 'http://www.w3.org/2001/xml-events',
        'xmlns:xlink': 'http://www.w3.org/1999/xlink',
        # "fill": "#FFFFFF",
    }

    return textwrap.dedent(r'''
        <?xml version="1.0" encoding="utf-8" ?>
        <svg {attrs:s}>
        {data:s}
        </svg>
    ''').format(
        attrs=' '.join(['{key:s}="{val:s}"'.format(
            key=key, val=props[key]) for key in props]),
        data=area.svg()
    ).strip().replace('stroke-width="2.0"', 'stroke-width="0.0"')


def write_img(f, c, hidden_layer=[]):
    pgds = os.path.join(f"{f}.gds")
    c.write_gds(pgds, with_metadata=True)
    if hidden_layer:
        if hidden_layer[0] is int:
            hidden_layer = [hidden_layer]
    psvg = os.path.join(f"{f}.svg")
    gdstk.read_gds(pgds).top_level()[0].write_svg(
        psvg,
        scaling=1,
        background="#000000",
        shape_style={k: {"fill": "none",
                                 "stroke": "none",
                                 "stroke-dasharray": "8,8"} for k in hidden_layer},
        pad=0)
    return io.open(psvg, "r").read()
    # renderPM.drawToFile(svg2rlg("tmp"), f, fmt='PNG')


def pic2gds(fileName, sizeOfTheCell, layerNum=1, isDither=False, scale=1):
    main(fileName, sizeOfTheCell, layerNum, isDither, scale)
    return "image.bmp", "image.gds"


def finish(c, name):
    c.add_label(name, position=c.bbox_np()[1])


def normal_from_orientation(orientation):
    return [cos(orientation/180*pi), sin(orientation/180*pi)]


def raster_slice(scene, dx, center, w, h, normal):
    # scene.show()
    g = sum(scene.geometry.values())
    s = g.section(plane_origin=center, plane_normal=normal)
    tangent = np.array([[-normal[1], normal[0], 0]]).T
    normal = np.array(normal)
    normal.shape = (3, 1)
    center = np.array(center)
    center.shape = (3, 1)
    tangent2 = np.cross(normal.T, tangent.T).T

    R = np.linalg.inv(np.concatenate((tangent, tangent2, normal), axis=1))
    to_2D = np.eye(4)
    to_2D[:3, :3] = R

    v = -tangent*w/2-h/2*tangent2
    xmin,    ymin, _ = numpy.reshape(np.matmul(R, v), (3,))

    v = +tangent*w/2+h/2*tangent2
    xmax,    ymax, _ = numpy.reshape(np.matmul(R, v), (3,))

    xmin, xmax = min(xmin, xmax), max(xmax, xmin)
    ymin, ymax = min(ymin, ymax), max(ymax, ymin)

    x, y, z = numpy.reshape(np.matmul(R, center), (3,))
    xmin, xmax = xmin+x, xmax+x
    ymin, ymax = ymin+y, ymax+y

    if not s:
        return np.zeros((round(h/dx), round(w/dx),)).T
    p, _ = s.to_planar(to_2D=to_2D)
    p = p.polygons_full

    # fig = plt.figure(1, )
    # ax = fig.add_subplot(121)
    # plot_polygon(p[0])
    # plt.show()

    p = clip_by_rect(GeometryCollection(p),
                     xmin, ymin, xmax, ymax)
    #  x-w/2, y-h/2, x + w/2, y + h/2)
    if not p:
        return None
    # plot_polygon(p,)
    # plt.show()

    p = shapely.affinity.translate(p, (xmax-xmin)/2-x, (ymax-ymin)/2-y)
    svg = s2svg(p, [0, 0, w, h], dx)
    svg2png(bytestring=svg, write_to='temp.png', background_color="#00000000",
            output_width=round(w/dx), output_height=round(h/dx))
    # p.plot()
    # with open('_.svg', 'w') as f:
    #     f.write(svg)
    # plot_polygon(p,)
    # plt.show()

    img = iio.imread('temp.png')
    img = np.sum(img, 2).T
    d = np.max(img)-np.min(img)
    b = np.min(img)
    if d == 0:
        if b == 0:
            d = 1
        else:
            d = b
            b = 0
    img = (img-b)/(d)
    # img = rasterio.features.rasterize(
    #     [p], transform=(dx, 0, 0, 0, dx, 0, 0, 0, dx), out_shape=(round(h/dx), round(w/dx),)).T
    # plt.imshow(img)
    # plt.show()

    return img


p = shapely.Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])
svg = s2svg(p, [0, 0, 2, 2], .2)
svg2png(bytestring=svg, write_to='temp.png', background_color="#00000000",
        output_width=10, output_height=10)


def material_slice(c, dx, center, w, h, normal, layers, layer_stack, materials=MATERIALS):
    LAYER_VIEWS["WGCLAD"].visible = True
    view = LAYER_VIEWS.layer_views["WGCLAD"]

    _layer_views = copy.deepcopy(LAYER_VIEWS)
    _layer_views.layer_views.clear()
    for i, layer in enumerate(layers):
        k = f"l{i}"
        _layer_views.layer_views[k] = copy.deepcopy(view)
        _layer_views.layer_views[k].layer = layer

    layers = sum([[[v.mesh_order, v.material, tuple(layer), k]
                 for k, v in get_layers(layer_stack, layer, withkey=True)] for layer in layers], [])

    epsmin = 100
    eps_array = 0
    layers = sorted(layers, key=lambda x: -x[0])
    for layer in layers:
        m = layer[1]
        k = layer[3]
        eps = materials[m]["epsilon"]

        _layer_stack = copy.deepcopy(layer_stack)
        _layer_stack.layers.clear()
        _layer_stack.layers[k] = layer_stack.layers[k]
        scene = c.to_3d(
            layer_stack=_layer_stack,
            layer_views=_layer_views)
        # mesh = sum(scene.geometry.values())
        _mask = raster_slice(scene,
                             dx, w=w, h=h,
                             center=center,
                             normal=normal)
        # print(layer, center)
        if _mask is None:
            # print("no mask ")
            pass
        else:
            if eps < epsmin:
                epsmin = eps
            eps_array = _mask*eps+(1-_mask)*eps_array
    holes = np.where(eps_array < epsmin, 1, 0)
    eps_array = (1-holes)*eps_array+epsmin*holes
    # plt.clf()
    # plt.imshow(eps_array)
    # plt.show()
    return eps_array


def material_voxelate(c, dl, zmin, zmax, layers, layer_stack, path):
    stacks = sum([[[v.mesh_order, v.material, tuple(layer), k]
                 for k, v in get_layers(layer_stack, layer, withkey=True)] for layer in layers], [])

    stacks = sorted(stacks, key=lambda x: -x[0])
    layer_stack_info = dict()
    for i, stack in enumerate(stacks):
        m = stack[1]
        l1, l2 = layer = stack[2]
        k = stack[3]
        # eps = materials[m]["epsilon"]

        _layer_stack = copy.deepcopy(layer_stack)
        _layer_stack.layers.clear()

        d = copy.deepcopy(layer_stack.layers[k])
        if d.zmin <= zmax and d.bounds[1] >= zmin:
            _layer_stack.layers[k] = d
            d.zmin = max(zmin, d.zmin)
            d.thickness = min(zmax-d.zmin, d.thickness)
            # _d.bounds = (_d.zmin, _d.zmin+_d.thickness)
            origin = (c.extract([layer]).bbox_np()-c.bbox_np())[0].tolist()
            gf.export.to_stl(c, os.path.join(
                path, f"{k}.stl"), layer_stack=_layer_stack)
            dir = os.path.join(path, k)
            os.makedirs(dir, exist_ok=True)
            stltovoxel.convert_file(
                os.path.join(path, f'{k}_{l1}_{l2}.stl'), os.path.join(dir, 'output.png'), voxel_size=dl, pad=0)
            layer_stack_info[k] = {
                "layer": (l1, l2),
                "zmin": d.zmin,
                "thickness": d.thickness,
                "material": m,
                "mesh_order": stack[0],
                "origin": origin,
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
