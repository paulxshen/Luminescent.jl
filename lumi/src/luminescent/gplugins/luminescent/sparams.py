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
                    dx=.05,
                    entries=[],
                    wavelength=1.55, center_wavelength=None, keys=[],
                    approx_2D=False, layer_stack=LAYER_STACK,
                    study="sparams",
                    **kwargs):

    d = 2 if approx_2D else 3
    prob = dict()
    ports = [int(p.name[1])
             for p in c.get_ports_list(prefix="o")]

    if not entries:
        if not keys:
            keys = []
            for i in ports:
                for o in ports:
                    keys.append(f"o,i")

        if type(wavelength) not in [list, tuple]:
            wavelengths = [wavelength]
        else:
            wavelengths = wavelength
        wavelengths = sorted(wavelengths)
        for w in wavelengths:
            for k in keys:
                entries.append([w, *unpack_sparam_key(k)])

    l = []
    for w, po, mo, pi, mi in entries:
        k = [w, pi, mi, pi, mi]
        if k not in entries:
            l.append(k)
    entries.extend(l)

    wio = {}
    for w, po, mo, pi, mi in entries:
        if w not in wio:
            wio[w] = {}

        if pi not in wio[w]:
            wio[w][pi] = {"mn": mi, "o": {}}
        else:
            wio[w][pi]["mn"] = max(wio[w][pi]["mn"], mi)

        if po not in wio[w][pi]["o"]:
            wio[w][pi]["o"][po] = mo
        else:
            wio[w][pi]["o"][po] = max(wio[w][pi]["o"][po], mo)

    runs = []
    for w in wio:
        runs.extend([{
            "d": d,
            "sources": {
                port_number(i): {
                    "center": (np.array(c.ports[f"o{port_number(i)}"].center)/1e3).tolist(),
                    "width": (np.array(c.ports[f"o{port_number(i)}"].width)/1e3).tolist(),
                    "normal": normal_from_orientation(c.ports[f"o{port_number(i)}"].orientation),
                    # "endpoints": extend(c.ports[f"o{port_number(i)}"]., margin),
                    "wavelength_mode_numbers": {w: range(wio[w][i]["mn"]+1) for w in wavelengths},
                    "port": port_number(i),
                }},
            "monitors": {
                port_number(o): {
                    "port": port_number(o),
                    "normal": normal_from_orientation(c.ports[f"o{port_number(o)}"].orientation),
                    "center": (np.array(c.ports[f"o{port_number(o)}"].center)/1e3).tolist(),
                    "width": (np.array(c.ports[f"o{port_number(o)}"].width)/1e3).tolist(),
                    # "endpoints": extend(c.ports[f"o{port_number(o)}"].endpoints, margin),
                    "wavelength_mode_numbers": {w: range(wio[w][i]["o"][o]+1) for w in wavelengths},
                } for o in wio[i]["o"]},
        } for i in wio[w]])

    prob = setup(c, study=study,  dx=dx,
                 runs=runs, margin=margin,
                 layer_stack=layer_stack, approx_2D=approx_2D, **kwargs)
    prob["wavelengths"] = wavelengths

    return prob

    # l = [k for k in wio[w] if port_number(k) == pi]
    # if not l:
    #     wio[w][f"o{pi}@{mi}"] = []
    # else:
    #     k = l[0]
    #     mn = max(mode_number(k), mi)
    #     if mn != mode_number(k):
    #         wio[w][i] = wio[w][k]
    #         del wio[w][k]

    # l = [k for k in wio[w][i] if port_number(k) == po]
    # if not l:
    #     wio[w][f"o{pi}@{mi}"]
    # else:
    #     k = l[0]
    #     mn = max(mode_number(k), mi)
    #     if mn != mode_number(k):
    #         wio[w][f"o{pi}@{mn}"] = wio[w][k]
    #         del wio[w][k]

    # if po not in wio[w][pi]:
    #     wio[w][pi]["o"][po] = mo
    # else:
    #     wio[w][pi]["o"][po] = max(wio[w][pi]["o"][po], mo)
