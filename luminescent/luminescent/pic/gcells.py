
from ..constants import *
from ..utils import *
import gdsfactory as gf
from gdsfactory.cross_section import Section
from networkx import center
from numpy import cumsum
from gdsfactory.generic_tech import LAYER_STACK, LAYER

# import .utils as utils
# from layers import *


def mimo(west=0, east=0, south=0, north=0,
         l=2.0, w=2.0, wwg=.5, lwg=None, taper=0,
         wwg_west=None, wwg_east=None, wwg_south=None, wwg_north=None,
         wwg_layer=LAYER.WG,  # bbox_layer=LAYER.WAFER,
         design_layer=DESIGN_LAYER,
         **kwargs):
    design = gf.Component()
    c = gf.Component(**kwargs)
    if lwg is None:
        lwg = 4*wwg
    p = [(0, 0), (l, 0), (l, w), (0, w)]
    design.add_polygon(p,                       layer=design_layer)
    c.add_polygon(p,                       layer=wwg_layer)

    ld = [west,  east, south, north]
    for i, v, d in zip(range(4), ld, [w, w, l, l]):
        if type(v) is int:
            ld[i] = [(.5+j)*d/v for j in range(v)]
    lwwg = [wwg_west, wwg_east, wwg_south, wwg_north]
    for i, v in enumerate(lwwg):
        if v is None:
            v = wwg
        if type(v) is float or type(v) is int:
            lwwg[i] = [v]*len(ld[i])

    n = 0
    for (i, x, y, d, wwg, a) in zip(
        range(4),
        [0,  l, 0, 0],
        [0, 0, 0, w],
        ld,
        lwwg,
        [180, 0, -90, 90]
    ):
        for wwg, v in zip(wwg, d):
            center = (x, y+v) if i in [0, 1] else (x+v, y)
            wwg2 = wwg+2*taper*lwg
            name = "o"+str(n+1)
            design.add_port(name, center=center, width=wwg2,
                            orientation=a, layer=wwg_layer)
            wg = c << gf.components.taper(
                length=lwg, width1=wwg, width2=wwg2, layer=wwg_layer)
            wg.connect("o2", design.ports[name])
            c.add_port(name, port=wg.ports["o1"])
            n += 1

    design = c << design
    return c
