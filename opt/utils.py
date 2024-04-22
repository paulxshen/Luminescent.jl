# from svglib.svglib import svg2rlg
# from reportlab.graphics import renderPM
# import cairosvg
import io
import gdstk
import pylab
import EMpy
# import EMpy as em
import numpy
from functools import partial
from math import cos, pi, sin

import matplotlib.pyplot as plt
import numpy as np

from gdsfactory.cross_section import Section
import gdsfactory as gf
import bson
epscore = 3.4757**2
epsclad = 1.444**2


def solve_modes(wwg, hwg, λ=1.55, wm=.1, hm=.1, dx=.05, epscore=epscore, epsclad=epsclad):
    _ = 1e-9

    w = wwg+2*wm
    h = hwg+2*hm
    x = numpy.linspace(0, w-dx, round(w/dx))
    y = numpy.linspace(0, h-dx, round(h/dx))

    def epsfunc(x_, y_):
        """Return a matrix describing a 2d material.

        :param x_: x values
        :param y_: y values
        :return: 2d-matrix
        """
        xx, yy = numpy.meshgrid(x_, y_)
        return numpy.where(
            # (numpy.abs(xx.T - 1.24) <= 0.24) * (numpy.abs(yy.T - 1.11) <= 0.11),
            (xx.T < w-wm-_) * (xx.T > wm-_)*(yy.T > hm-_)*(yy.T < h-hm-_),
            epscore,
            epsclad,
        )

    neigs = 1
    tol = 1e-8
    boundary = "0000"

    eps=epsfunc(x, y)
    solver = EMpy.modesolvers.FD.VFDModeSolver(λ, x, y, epsfunc, boundary).solve(
        neigs, tol
    )

    fig = pylab.figure()

    fig.add_subplot(1, 3, 1)
    Ex = numpy.transpose(solver.modes[0].get_field("Ex", x, y))
    pylab.contourf(x, y, abs(Ex), 50)
    pylab.title("Ex")

    modes = [{k: m.get_field(k, x, y) for k in [
        "Ex", "Ey", "Ez", "Hx", "Hy", "Hz"]} for m in solver.modes]
    # mode = {"Ex": Ex, "Ey": Ey, "Ez": Ez, "Hx": Hx, "Hy": Hy, "Hz": Hz}
    # np.savez("modes_{wwg}_{hwg}.npz", **mode)
    pylab.show()
    # return solver.modes
    return modes,eps


def write_img(f, c, hidden_layer=[]):
    c.write_gds(f"{f}.gds", ".")
    if hidden_layer:
        if hidden_layer[0] is int:
            hidden_layer = [hidden_layer]
    gdstk.read_gds("tmp").top_level()[0].write_svg(
        f"{f}.svg",
        shape_style={k: {"fill": "none",
                                 "stroke": "none",
                                 "stroke-dasharray": "8,8"} for k in hidden_layer},
        pad=0)
    # png = cairosvg.svg2png(url="tmp", write_to=f)

    # read svg -> write png
    # renderPM.drawToFile(svg2rlg("tmp"), f, fmt='PNG')
