
import cv2
import rasterio.features
import shapely
from shapely import GeometryCollection, Polygon
from shapely.ops import clip_by_rect
from shapely.plotting import plot_polygon, plot_points
import picToGDS
import io
import os
import sys
import gdstk
import pylab
import EMpy
# import EMpy as em
import numpy
from functools import partial
from math import cos, pi, sin, tan
import trimesh
import matplotlib.pyplot as plt
import numpy as np

from gdsfactory.cross_section import Section
import gdsfactory as gf
import bson

material_eps = {
    "si": 3.4757**2,
    "sio2": 1.444**2,
}


def solve_modes(eps, λ, dx):
    tol = 1e-4
    m, n = eps.shape
    x = numpy.linspace(0, (m-1)*dx, m)
    y = numpy.linspace(0, (n-1)*dx, n)

    def ϵfunc(x_, y_):
        return cv2.resize(eps, dsize=(len(x_), len(y_)), interpolation=cv2.INTER_CUBIC).T

    neigs = 1
    boundary = "0000"

    tol = 1e-6
    solver = EMpy.modesolvers.FD.VFDModeSolver(λ, x, y, ϵfunc, boundary).solve(
        neigs, tol
    )

    fig = pylab.figure()

    fig.add_subplot(1, 3, 1)
    Ex = numpy.transpose(solver.modes[0].get_field("Ex", x, y))
    pylab.imshow(abs(Ex))
    pylab.title("Ex")

    fig.add_subplot(1, 3, 2)
    Hy = numpy.transpose(solver.modes[0].get_field("Hy", x, y))
    pylab.imshow(abs(Hy))
    pylab.title("Hy")
    fig.add_subplot(1, 3, 3)
    e = numpy.transpose(ϵfunc(x, y))
    pylab.imshow(e)
    pylab.title("eps")

    modes = [{k: m.get_field(k, x, y) for k in [
        "Ex", "Ey", "Ez", "Hx", "Hy", "Hz"]} for m in solver.modes]
    return modes

# def solve_modes(wwg, hwg, λ=1.55, wm=.15, hm=.15, dx=.05, ϵcore=ϵcore, ϵclad=ϵclad):
    # tol = 1e-4

    # w = wwg+2*wm
    # h = hwg+2*hm
    # # x = numpy.linspace(0, w, round(w/dx)+1)
    # # y = numpy.linspace(0, h, round(h/dx)+1)
    # x = numpy.linspace(0, w, round(w/dx)+1)
    # y = numpy.linspace(0, h, round(h/dx)+1)

    # def ϵfunc(x_, y_):
    #     """Return a matrix describing a 2d material.

    #     :param x_: x values
    #     :param y_: y values
    #     :return: 2d-matrix
    #     """
    #     # xx, yy = numpy.meshgrid(x_, y_)
    #     # return numpy.where(
    #     #     # (numpy.abs(xx.T - 1.24) <= 0.24) * (numpy.abs(yy.T - 1.11) <= 0.11),
    #     #     (xx.T < w-wm-tol) * (xx.T > wm-tol) * \
    #     #     (yy.T > hm-tol)*(yy.T < h-hm-tol),
    #     #     ϵcore,
    #     #     ϵclad,
    #     # )

    #     def f(x, y):
    #         if (x < w-wm-tol) and (x > wm+tol) and (y > hm+tol) and (y < h-hm-tol):
    #             return ϵcore
    #         elif (x > w-wm+tol) or (x < wm-tol) or (y < hm-tol) or (y > h-hm+tol):
    #             return ϵclad
    #         else:
    #             return (ϵcore+ϵclad)/2
    #     # return np.ma[f(x, y) for x, y in zip(xx, yy)]
    #     r = np.array([[f(x, y) for x in x_] for y in y_]).T
    #     # r = np.array([1 for x in 1:2 for y in 1:3])
    #     return r
    # neigs = 1
    # tol = 1e-8
    # boundary = "0000"

    # ϵ = ϵfunc(x, y)
    # solver = EMpy.modesolvers.FD.VFDModeSolver(λ, x, y, ϵfunc, boundary).solve(
    #     neigs, tol
    # )

    # fig = pylab.figure()

    # fig.add_subplot(1, 3, 1)
    # Ex = numpy.transpose(solver.modes[0].get_field("Ex", x, y))
    # pylab.imshow(abs(Ex))
    # pylab.title("Ex")

    # fig.add_subplot(1, 3, 2)
    # Hy = numpy.transpose(solver.modes[0].get_field("Hy", x, y))
    # pylab.imshow(abs(Hy))
    # pylab.title("Hy")
    # fig.add_subplot(1, 3, 3)
    # e = numpy.transpose(ϵfunc(x, y))
    # pylab.imshow(e)
    # pylab.title("eps")

    # modes = [{k: m.get_field(k, x, y) for k in [
    #     "Ex", "Ey", "Ez", "Hx", "Hy", "Hz"]} for m in solver.modes]
    # # mode = {"Ex": Ex, "Ey": Ey, "Ez": Ez, "Hx": Hx, "Hy": Hy, "Hz": Hz}
    # # np.savez("modes_{wwg}tol{hwg}.npz", **mode)
    # # pylab.show()
    # # return solver.modes
    # return modes, ϵ, [w, h]


def write_img(f, c, hidden_layer=[]):
    pgds = os.path.join(f"{f}.gds")
    c.write_gds(pgds, with_metadata=True)
    if hidden_layer:
        if hidden_layer[0] is int:
            hidden_layer = [hidden_layer]
    psvg = os.path.join(f"{f}.svg")
    gdstk.read_gds(pgds).top_level()[0].write_svg(
        psvg,
        scaling=1000,
        background="#000000",
        shape_style={k: {"fill": "none",
                                 "stroke": "none",
                                 "stroke-dasharray": "8,8"} for k in hidden_layer},
        pad=0)
    return io.open(psvg, "r").read()
    # png = cairosvg.svg2png(url="tmp", write_to=f)

    # read svg -> write png
    # renderPM.drawToFile(svg2rlg("tmp"), f, fmt='PNG')


def pic2gds(fileName, sizeOfTheCell, layerNum=1, isDither=False, scale=1):
    picToGDS.main(fileName, sizeOfTheCell, layerNum, isDither, scale)
    return "image.bmp", "image.gds"


def finish(c, name):
    c.add_label(name, position=c.bbox[1])


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

    p, _ = s.to_planar(to_2D=to_2D)
    p = p.polygons_full
    if not p:
        return None

    # fig = plt.figure(1, )
    # ax = fig.add_subplot(121)
    plot_polygon(p[0])
    # plt.show()

    p = clip_by_rect(GeometryCollection(p),
                     xmin, ymin, xmax, ymax)
    #  x-w/2, y-h/2, x + w/2, y + h/2)
    if not p:
        return None
    plot_polygon(p,)
    # plt.show()

    p = shapely.affinity.translate(p, (xmax-xmin)/2-x, (ymax-ymin)/2-y)
    plot_polygon(p,)
    # plt.show()

    img = rasterio.features.rasterize(
        [p], transform=(dx, 0, 0, 0, dx, 0, 0, 0, dx), out_shape=(round(h/dx)+1, round(w/dx)+1,))
    plt.imshow(img)
    # plt.show()

    return img
