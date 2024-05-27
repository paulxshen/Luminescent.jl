# from gcells import ubend
from operator import inv
from cv2 import solve
from regex import E
from gcells import *
from prob import *
import pcells
import gdsfactory as gf
# import picToGDS as pg
import utils
from gdsfactory.generic_tech import LAYER, LAYER_STACK, LayerMap
from gdsfactory.technology import (
    LayerLevel,
    LayerStack,
    LayerView,
    LayerViews,
    LayerMap,
)
from gdsfactory.typings import Layer

dx = .05
wwg = 0.4
r = .4
h = 0.2
l = 3
lmin = 0.1


class MyLayerMap(LayerMap):
    WG: Layer = (1, 0)
    WGCLAD: Layer = (111, 0)
    BOX: Layer = (112, 0)
    PORT: Layer = (1, 10)


LAYER = MyLayerMap()

box_thickness = 0.5
thickness_clad = 0.5
LAYER_STACK.layers.update(dict(
    box=LayerLevel(
        layer=LAYER.BOX,
        thickness=box_thickness,
        zmin=-box_thickness,
        material="sio2",
        mesh_order=9,
    ),
    clad=LayerLevel(
        layer=LAYER.WGCLAD,
        zmin=0.0,
        material="sio2",
        thickness=thickness_clad,
        mesh_order=10,
    ),))
c, signals, targets, design = ubend_template(
    r,  l, wwg, dx=dx, wavelength=1.55, layer_map=LAYER)
prob = inverse_design_problem(
    c, lmin=lmin, dx=dx, signals=signals, targets=targets, design=design, layer_stack=LAYER_STACK)
solve(prob)
raise Exception("Not implemented")
utils.pic2gds("opt/ubend.png", dx)
init = gf.import_gds("image.gds", "TOP",  read_metadata=True).rotate(90)
ubend = ubend_template(r, wwg, h, l, dx=dx, dir="opt", lmin=lmin, init=init)
ubend.show()
ps = pcells.phase_shifter(l=20, wwg=wwg, nbends=4, ubend=ubend)
res = pcells.ubend_resonator(l=50, wwg=wwg, ubend=ubend)
ps.show()
res.show()
