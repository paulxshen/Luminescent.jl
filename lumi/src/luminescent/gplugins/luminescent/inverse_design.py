from .materials import MATERIALS
from .sparams import *
from .setup import *
from .constants import *
from .layers import *
from .utils import *
import bson
import gdsfactory as gf
from copy import deepcopy
from functools import partial
from math import cos, pi, sin
import os
from gdsfactory.generic_tech import LAYER_STACK, LAYER


def gcell_problem(c,  targets, iters,
                  lvoid=0, lsolid=0, symmetries=[],
                  weights=dict(),
                  eta=.1, init=1,   stoploss=.03,
                  design_region_layer=DESIGN_LAYER,
                  #    design_guess_layer=LAYER.GUESS,
                  fill_layer=LAYER.WG,
                  void_layer=None,
                  layer_stack=LAYER_STACK, materials=MATERIALS,
                  plot=False, approx_2D=True,
                  restart=True, save_memory=False, **kwargs):
    design_region_layer = tuple(design_region_layer)
    # if not approx_2D:
    #     raise NotImplementedError(
    #         "3D inverse design feature must be requested from Luminescent AI info@luminescentai.com")

    eta *= 10
    if "phase_shifter" in targets:
        keys = ["o2@0,o1@0"]
        wavelengths = [preset["wavelength"]]
    else:
        if "tparams" in targets:
            for wl in targets["tparams"]:
                d = {}
                for k in targets["tparams"][wl]:
                    po, mo, pi, mi = unpack_sparam_key(k)
                    if mo == 0:
                        _mo = 1
                    elif mo == 1:
                        _mo = 0
                    else:
                        mo = None

                    if mo is not None:
                        _k = f"o{po}@{_mo},o{pi}@{mi}"
                        if _k not in targets["tparams"][wl]:
                            d[_k] = 0

                    _k = f"o{pi}@0,o{pi}@{mi}"
                    if _k not in targets["tparams"][wl]:
                        d[_k] = 0
                    _k = f"o{pi}@1,o{pi}@{mi}"
                    if _k not in targets["tparams"][wl]:
                        d[_k] = 0
                targets["tparams"][wl].update(d)
        print(targets)

        targets1 = {}
        keys = SortedSet()
        wavelengths = SortedSet()
        for k, d in targets.items():
            d = {
                wl: {
                    long_sparam_key(k): v for k, v in d.items()
                } for wl, d in d.items()}
            targets1[k] = d
            if k in ["sparams", "tparams"]:
                for s in sum([list(d.keys()) for d in d.values()], []):
                    keys.add(s)
            for wl in d:
                wavelengths.add(wl)
        targets = targets1
        keys = list(keys)
        print(keys)
        wavelengths = list(wavelengths)

    prob = sparams_problem(c,
                           layer_stack=layer_stack, materials=materials,
                           wavelength=wavelengths,
                           study="inverse_design",
                           keys=keys,
                           approx_2D=approx_2D, ** kwargs)
    prob["restart"] = restart
    prob["weights"] = {**{
        "tparams": 1,
        "sparams": 1,
        "phasediff": 1,
    }, **weights}
    prob["save_memory"] = save_memory
    prob["targets"] = targets
    prob["wavelengths"] = wavelengths
    # prob["init"] = init
    prob = {**prob, **kwargs}
    prob["eta"] = eta
    polys = c.extract([design_region_layer]).get_polygons()
    if not symmetries:
        symmetries = [[]]*len(polys)
    else:
        if type(symmetries[0]) is not list:
            symmetries = [symmetries]*len(polys)

    def _bbox(b):
        return [[b.left/1e3, b.bottom/1e3], [b.right/1e3, b.top/1e3]]
    prob["designs"] = [
        {
            "layer": design_region_layer,
            "bbox": _bbox(p.bbox()),
            "symmetries": s,
            "init": init,
            "lvoid": lvoid,
            "lsolid": lsolid,
        } for p, s in zip(list(polys.values())[0], symmetries)
    ]
    epsmin = np.min(prob["eps_2D"])
    prob["stoploss"] = stoploss
    prob["design_config"] = dict()
    l = get_layers(layer_stack, fill_layer)[0]
    d = {"thickness": l.thickness, "material": l.material, "zmin": l.zmin}
    d["epsilon"] = materials[d["material"]]["epsilon"]
    d["ϵ"] = d["epsilon"]
    d["layer"] = fill_layer
    prob["design_config"]["fill"] = d

    if void_layer is not None:
        l = get_layers(layer_stack, void_layer)[0]
        d = {"thickness": l.thickness, "material": l.material, "zmin": l.zmin}
        d["layer"] = void_layer
        d["epsilon"] = materials[d["material"]]["epsilon"]
        d["ϵ"] = d["epsilon"]
    else:
        d = copy.deepcopy(d)
        d["epsilon"] = epsmin
        d["ϵ"] = d["epsilon"]
    prob["design_config"]["void"] = d

    prob["design_config"]["design_region_layer"] = design_region_layer
    prob["iters"] = iters
    return prob


def apply_design(c0,  sol):
    path = sol["path"]
    a = gf.Component()
    a.add_ref(c0)
    fill = sol["design_config"]["fill"]["layer"]
    dl = sol["design_config"]["design_region_layer"]
    for i, d in enumerate(sol["designs"]):
        x0, y0 = d["bbox"][0]
        x1, y1 = d["bbox"][1]
        b = gf.Component()
        b.add_polygon(
            [(x0, y0), (x1, y0), (x1, y1), (x0, y1)], layer=fill)
        a = gf.boolean(a, b, "not", layer=fill)
    c = gf.Component()
    c << a
    # for layer in c0.layers:
    #     if layer != dl:
    #         c.add_ref(c0.extract([layer]))
    # c.show()
    # raise ValueError()
    g = gf.import_gds(os.path.join(path, f"optimized_design_region_{i+1}.gds"))
    polygons = g.get_polygons(merge=True)
    g = gf.Component()
    for p in polygons[1]:
        g.add_polygon(p, layer=fill)
    g = c << g
    g.xmin = x0
    g.ymin = y0
    return c
