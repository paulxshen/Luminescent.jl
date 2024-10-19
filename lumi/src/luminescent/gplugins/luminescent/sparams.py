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
    ports = [p.name for p in c.get_ports_list(prefix="o")]

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
    wavelengths = SortedSet()
    for w, po, mo, pi, mi in entries:
        k = [w, pi, mi, pi, mi]
        if k not in entries:
            l.append(k)
        wavelengths.add(w)
    entries.extend(l)
    wavelengths = list(wavelengths)

    wimo = {}
    for w, po, mo, pi, mi in entries:
        if w not in wimo:
            wimo[w] = {}
        if pi not in wimo[w]:
            wimo[w][pi] = {}
        if mi not in wimo[w][pi]:
            wimo[w][pi][mi] = {}

        if po not in wimo[w][pi][mi]:
            wimo[w][pi][mi][po] = mo
        else:
            wimo[w][pi][mi][po] = max(wimo[w][pi][mi][po], mo)

    runs = []
    for w in wimo:
        for i in wimo[w]:
            for mi in wimo[w][i]:
                runs.append({
                    "sources": {
                        i: {
                            "center": (np.array(c.ports[i].center)/1e3).tolist(),
                            "width": (np.array(c.ports[i].width)/1e3).tolist(),
                            "normal": normal_from_orientation(c.ports[i].orientation),
                            "wavelength_mode_numbers": {w: [mi]},
                            "port": i,
                        }},
                    "monitors": {
                        o: {
                            "port": o,
                            "normal": normal_from_orientation(c.ports[o].orientation),
                            "center": (np.array(c.ports[o].center)/1e3).tolist(),
                            "width": (np.array(c.ports[o].width)/1e3).tolist(),
                            # "endpoints": extend(c.ports[o].endpoints, margin),
                            "wavelength_mode_numbers": {w: list(range(wimo[w][i][mi][o]+1))}
                        } for o in wimo[w][i][mi]}})

    prob = setup(c, study=study,  dx=dx,
                 runs=runs, margin=margin,
                 layer_stack=layer_stack, approx_2D=approx_2D, **kwargs)
    prob["wavelengths"] = wavelengths

    return prob

    # l = [k for k in wimo[w] if port_number(k) == pi]
    # if not l:
    #     wimo[w][f"o{pi}@{mi}"] = []
    # else:
    #     k = l[0]
    #     mn = max(mode_number(k), mi)
    #     if mn != mode_number(k):
    #         wimo[w][i] = wimo[w][k]
    #         del wimo[w][k]

    # l = [k for k in wimo[w][i] if port_number(k) == po]
    # if not l:
    #     wimo[w][f"o{pi}@{mi}"]
    # else:
    #     k = l[0]
    #     mn = max(mode_number(k), mi)
    #     if mn != mode_number(k):
    #         wimo[w][f"o{pi}@{mn}"] = wimo[w][k]
    #         del wimo[w][k]

    # if po not in wimo[w][pi]:
    #     wimo[w][pi]["o"][po] = mo
    # else:
    #     wimo[w][pi]["o"][po] = max(wimo[w][pi]["o"][po], mo)
