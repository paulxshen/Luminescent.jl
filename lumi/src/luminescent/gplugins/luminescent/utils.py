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
from .generic_tech import LAYER_STACK, LAYER, LAYER_VIEWS
from .materials import MATERIAL_LIBRARY
from .constants import PATH
from gdsfactory.cross_section import Section
import gdsfactory as gf
import bson


def trim(x, dx):
    return round(x/dx)*dx


def extend(endpoints, wm):
    v = endpoints[1]-endpoints[0]
    v = v/np.linalg.norm(v)
    return [
        (endpoints[0]-wm*v).tolist(), (endpoints[1]+wm*v).tolist()]


def portsides(c):
    res = [[False, False], [False, False]]
    bbox = c.bbox_np()
    xmin0, ymin0 = bbox[0]
    xmax0, ymax0 = bbox[1]
    for p in c.ports:
        x, y = np.array(p.center)/1e3
        if x == xmin0:
            res[0][0] = True
        if x == xmax0:
            res[1][0] = True
        if y == ymin0:
            res[0][1] = True
        if y == ymax0:
            res[1][1] = True
    return res


def add_bbox(c, layers, nonport_margin=0):
    margin = nonport_margin
    bbox = c.bbox_np()
    xmin0, ymin0 = bbox[0]
    xmax0, ymax0 = bbox[1]
    xmin, ymin, xmax, ymax = xmin0-margin, ymin0-margin, xmax0+margin, ymax0+margin
    for p in c.ports:
        # p = c.ports[k]
        x, y = np.array(p.center)/1e3
        if x == xmin0:
            xmin = x
        if x == xmax0:
            xmax = x
        if y == ymin0:
            ymin = y
        if y == ymax0:
            ymax = y
    p = [(xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax)]
    _c = gf.Component()
    _c << c
    for layer in layers:
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
    solvera = EMpy.modesolvers.FD.VFDModeSolver(
        λ, x, y1, ϵfunc1,  "AA00").solve(neigs, tol)
    solvers = EMpy.modesolvers.FD.VFDModeSolver(
        λ, x, y1, ϵfunc1, "SS00").solve(neigs, tol)

    # if plot:
    # fig = pylab.figure()
    # fig.add_subplot(2, 3, 1)
    # Hy = numpy.transpose(solver.modes[0].get_field("Hy", x, y))
    # pylab.imshow(abs(Hy))
    # fig.add_subplot(2, 3, 2)
    # Hy = numpy.transpose(solvera.modes[0].get_field("Hy", x, y1))
    # pylab.imshow(abs(Hy))
    # fig.add_subplot(2, 3, 3)
    # Hy = numpy.transpose(solvers.modes[0].get_field("Hy", x, y1))
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

    v = sorted(solvera.modes, key=lambda x: -np.abs(x.neff))
    modes1 = [{k: m.get_field(k, x, y1)[:, 0] for k in [
        "Ex", "Ey", "Ez", "Hx", "Hy", "Hz"]} for m in v]
    neffs1 = [np.real(m.neff) for m in v]
    # "Ex", "Ey", "Ez", "Hx", "Hy", "Hz"]} for m in sorted(solvera.modes+solvers.modes, key=lambda x: -np.abs(x.neff))]
    # print(solver.modes[0].get_field("Ex", x, y).shape)
    return modes, modes1, neffs, neffs1


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
    svg2png(bytestring=svg, write_to='_.png', background_color="#00000000",
            output_width=round(w/dx), output_height=round(h/dx))
    # p.plot()
    # with open('_.svg', 'w') as f:
    #     f.write(svg)
    # plot_polygon(p,)
    # plt.show()

    img = iio.imread('_.png')
    img = np.sum(img, 2).T
    img = (img-np.min(img))/(np.max(img)-np.min(img))
    # img = rasterio.features.rasterize(
    #     [p], transform=(dx, 0, 0, 0, dx, 0, 0, 0, dx), out_shape=(round(h/dx), round(w/dx),)).T
    # plt.imshow(img)
    # plt.show()

    return img


p = shapely.Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])
svg = s2svg(p, [0, 0, 2, 2], .2)
svg2png(bytestring=svg, write_to='_.png', background_color="#00000000",
        output_width=10, output_height=10)


def material_slice(c, dx, center, w, h, normal, layers, layer_stack, layer_views=LAYER_VIEWS):
    layer_views = layer_views.copy()
    layer_views.layer_views["WGCLAD"].visible = True
    epsmin = 100
    eps_array = 0
    for layer in layers:
        layer_views.layer_views[f"{layer}"] = layer_views.layer_views["WGCLAD"]
        layer_views.layer_views[f"{layer}"].layer = layer

        # print(tuple(layer))
        m = get_layer(layer_stack, layer).material
        if m is not None:
            eps = MATERIAL_LIBRARY[m].epsilon
            scene = c.to_3d(
                layer_stack=layer_stack,
                layer_views=layer_views,
                exclude_layers=list(set(c.layers)-set([layer])))
            # mesh = sum(scene.geometry.values())
            _mask = raster_slice(scene,
                                 dx, w=w, h=h,
                                 center=center,
                                 normal=normal)
            # _mask = trimesh.voxel.creation.voxelize(mesh, dx)
            if eps < epsmin:
                epsmin = eps
            eps_array = _mask*eps+(1-_mask)*eps_array
    eps_array = eps_array+np.where(eps_array == 0, 1, 0)*epsmin
    # plt.clf()
    # plt.imshow(eps_array)
    # plt.show()
    return eps_array


def material_voxelate(c, dx, center, l, w, h,  normal, layers, layer_stack, layer_views=LAYER_VIEWS):
    return np.stack([material_slice(c, dx, [center[0]+x]+center[1:], w, h, normal,    layers, layer_stack, layer_views) for x in np.arange(dx/2, l-.001, dx)],)


def get_layer(layer_stack, layer):
    for x in layer_stack.layers.values():
        if hasattr(x.layer, "layer"):
            if tuple(x.layer.layer) == tuple(layer):
                return x

        # if hasattr(x.layer, "layer1") and x.layer.layer1.layer == layer:
        #     return x


def show_solution(path=None):

    if path is None:
        path = sorted(os.listdir(PATH))[-1]
        path = os.path.join(PATH, path)
    print(f"showing solution from {path}")
    # Load an image
    for s in ["after.png", "run1.png"]:
        try:
            img = Image.open(os.path.join(path, s))
        except:
            pass
        else:
            img.show()
            try:
                display(img)
            except:
                pass
            return
    print("no plot image found")
