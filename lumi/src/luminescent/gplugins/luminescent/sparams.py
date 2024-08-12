from .setup import *
from .constants import *
from .layers import *
from .utils import *
import gdsfactory as gf
from copy import deepcopy
import matplotlib.pyplot as plt
import numpy as np

from gdsfactory.cross_section import Section
from gdsfactory.generic_tech import LAYER_STACK, LAYER


def sparams_problem(c: gf.Component,
                    margin=None,  # zmargin=None,zlims=None,
                    dx=.05, wavelengths=[1.55], center_wavelength=None, keys=None,
                    approx_2D=False, layer_stack=LAYER_STACK,
                    study="sparams",
                    **kwargs):
    wavelengths = sorted(wavelengths)
    d = 2 if approx_2D else 3
    prob = dict()

    ports = [int(p.name[1])
             for p in c.get_ports_list(prefix="o")]
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

    if margin is None:
        margin = trim(2*max([p.width/1e3 for p in c.ports]), dx)
    runs = [{
        "d": d,
        "sources": {
            port_number(i): {
                "center": (np.array(c.ports[f"o{port_number(i)}"].center)/1e3).tolist(),
                "width": (np.array(c.ports[f"o{port_number(i)}"].width)/1e3+2*margin).tolist(),
                "normal": normal_from_orientation(c.ports[f"o{port_number(i)}"].orientation),
                # "endpoints": extend(c.ports[f"o{port_number(i)}"]., margin),
                "wavelength_mode_numbers": {w: [mode_number(i)] for w in wavelengths},
                "port": port_number(i),
            }},
        "monitors": {
            port_number(o): {
                "port": port_number(o),
                "normal": normal_from_orientation(c.ports[f"o{port_number(o)}"].orientation),
                "center": (np.array(c.ports[f"o{port_number(o)}"].center)/1e3).tolist(),
                "width": (np.array(c.ports[f"o{port_number(o)}"].width)/1e3+2*margin).tolist(),
                # "endpoints": extend(c.ports[f"o{port_number(o)}"].endpoints, margin),
                "wavelength_mode_numbers": {w: [mode_number(o)] for w in wavelengths},
            } for o in io[i]},
    } for i in io]

    prob = setup(c, study=study,  dx=dx,
                 runs=runs, margin=margin,
                 layer_stack=layer_stack, approx_2D=approx_2D, **kwargs)
    prob["wavelengths"] = wavelengths

    return prob
