from statistics import mode
from .setup import *
from ..constants import *
from ..layers import *
from ..utils import *
from ..materials import *
import gdsfactory as gf
from copy import deepcopy
import matplotlib.pyplot as plt
import numpy as np

from gdsfactory.cross_section import Section
from gdsfactory.generic_tech import LAYER_STACK, LAYER


def make_sim_prob(
    path,
    wavelengths=[],
    center_wavelength=None,
    frequencies=[],
    center_frequency=None,
        sources=[],
        monitors=[],
        nres=30,
        layer_stack={},
        materials={},
        study="simulation",
        dtype="float32",
        Ttrans=None,
        wl_res=.01,
        # Tss=None,
        gpu=None,
        margins=[[0, 0, 0], [0, 0, 0]],
):

    materials = {**MATERIALS, **materials}

    if not wavelengths:
        if not center_frequency:
            center_frequency = mode(frequencies)

        wavelengths = [center_wavelength * center_frequency /
                       f for f in sorted(frequencies, reverse=True)]
    wavelengths, center_wavelength, T = adjust_wavelengths(
        wavelengths, wl_res)

    dx = center_wavelength/nres
    dl = dx/4

    GEOMETRY = os.path.join(path, "geometry")
    bbox = [[None, None, None], [None, None, None]]
    for fn in os.listdir(GEOMETRY):
        if fn.lower().endswith(".stl"):
            STL = os.path.join(GEOMETRY, fn)
            pymeshfix.clean_from_file(STL, STL)
            mesh = pv.read(STL)
            for (i, v, w) in zip(range(3), bbox[0], mesh.bounds[0:3]):
                if v is None or w < v:
                    bbox[0][i] = w
            for (i, v, w) in zip(range(3), bbox[1], mesh.bounds[3:6]):
                if v is None or w > v:
                    bbox[1][i] = w
    bbox = [(np.array(b)-np.array(m)).tolist() for b, m in zip(bbox, margins)]

    for fn in os.listdir(GEOMETRY):
        if fn.lower().endswith(".stl"):
            STL = os.path.join(GEOMETRY, fn)
            mesh = pv.read(STL)
            im = stl_to_array(mesh, dl, bbox)
            name = fn[:-4]
            np.save(os.path.join(GEOMETRY, f'{name}.npy'), im)

    L = (bbox[1][0]-bbox[0][0], bbox[1][1]-bbox[0][1], bbox[1][2]-bbox[0][2])
    xs = np.linspace(bbox[0][0], bbox[1][0], 1 +
                     round((bbox[1][0]-bbox[0][0])/dl)).tolist()
    ys = np.linspace(bbox[0][1], bbox[1][1], 1 +
                     round((bbox[1][1]-bbox[0][1])/dl)).tolist()
    zs = np.linspace(bbox[0][2], bbox[1][2], 1 +
                     round((bbox[1][2]-bbox[0][2])/dl)).tolist()
#  dtype, center_wavelength, dl, xs, ys, zs, study, layer_stack, sources, monitors, materials, L, Ttrans, Tss, wavelengths
    prob = {
        "sources": sources,
        "monitors": monitors,
        "nres": nres,
        "study": study,
        'dtype': dtype,
        'center_wavelength': center_wavelength,
        'dl': dl,
        'dx': dx,
        'xs': xs,
        'ys': ys,
        'zs': zs,
        'layer_stack': layer_stack,
        'materials': materials,
        'L': L,
        'gpu_backend': gpu,
    }
    prob["wavelengths"] = wavelengths
    prob["Ttrans"] = None
    prob["Tss"] = T if len(wavelengths) > 1 else None
    PROB = os.path.join(path, "problem.json")
    with open(PROB, 'w') as f:
        json.dump(prob, f)
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
