from .setup import *
from ..constants import *
from ..layers import *
from ..utils import *
import gdsfactory as gf
from copy import deepcopy
import matplotlib.pyplot as plt
import numpy as np

from gdsfactory.cross_section import Section
from gdsfactory.generic_tech import LAYER_STACK, LAYER


def sparams_problem(c: gf.Component, nres,
                    wavelengths,
                    entries=None, keys=None,
                    N=3, layer_stack=LAYER_STACK,
                    study="sparams",
                    **kwargs):
    ports = [p.name for p in c.get_ports_list(prefix="o")]

    if not entries:
        entries = []
        if not keys:
            keys = []
            for i in ports:
                for o in ports:
                    keys.append(f"o,i")

        if type(wavelengths) not in [list, tuple, np.ndarray]:
            wavelengths = [wavelengths]
        else:
            wavelengths = wavelengths
        wavelengths, wl, T = adjust_wavelengths(wavelengths)
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

    imow = {}
    for w, po, mo, pi, mi in entries:
        if pi not in imow:
            imow[pi] = {}
        if mi not in imow[pi]:
            imow[pi][mi] = {}

        if po not in imow[pi][mi]:
            imow[pi][mi][po] = mo
        else:
            imow[pi][mi][po] = max(imow[pi][mi][po], mo)

    runs = []
    for _w in [1]:
        for i in imow:
            for mi in imow[i]:
                d = {
                    "sources": {
                        i: {
                            "center": (np.array(c.ports[i].center)/1e3).tolist(),
                            "width": (np.array(c.ports[i].width)/1e3).tolist(),
                            "normal": normal_from_orientation(c.ports[i].orientation),
                            "wavelength_mode_numbers": {w: [mi] for w in wavelengths},
                            "port": i,
                        }},
                    "monitors": {
                        o: {
                            "port": o,
                            "normal": normal_from_orientation(c.ports[o].orientation),
                            "center": (np.array(c.ports[o].center)/1e3).tolist(),
                            "width": (np.array(c.ports[o].width)/1e3).tolist(),
                            "wavelength_mode_numbers": {w: list(range(imow[i][mi][o]+1)) for w in wavelengths},
                        } for o in imow[i][mi]}}
                d["sources"] = SortedDict(d["sources"])
                d["monitors"] = SortedDict(d["monitors"])
                runs.append(d)

    prob = setup(c, study=study,  nres=nres, wl=wl,
                 runs=runs,
                 layer_stack=layer_stack, N=N, **kwargs)
    prob["wavelengths"] = wavelengths
    prob["Ttrans"] = None
    prob["Tss"] = T if len(wavelengths) > 1 else None
    return prob

    # l = [k for k in imow if port_number(k) == pi]
    # if not l:
    #     imow[f"o{pi}@{mi}"] = []
    # else:
    #     k = l[0]
    #     mn = max(mode_number(k), mi)
    #     if mn != mode_number(k):
    #         imow[i] = imow[k]
    #         del imow[k]

    # l = [k for k in imow[i] if port_number(k) == po]
    # if not l:
    #     imow[f"o{pi}@{mi}"]
    # else:
    #     k = l[0]
    #     mn = max(mode_number(k), mi)
    #     if mn != mode_number(k):
    #         imow[f"o{pi}@{mn}"] = imow[k]
    #         del imow[k]

    # if po not in imow[pi]:
    #     imow[pi]["o"][po] = mo
    # else:
    #     imow[pi]["o"][po] = max(imow[pi]["o"][po], mo)
