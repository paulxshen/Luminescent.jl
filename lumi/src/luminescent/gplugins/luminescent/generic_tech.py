import gdsfactory as gf
from gdsfactory.generic_tech.layer_map import LAYER
from gdsfactory.technology import LogicalLayer, LayerStack, LayerLevel
import gdsfactory as gf
from gdsfactory.generic_tech import LAYER, LAYER_STACK
from gdsfactory.generic_tech.get_klayout_pyxs import get_klayout_pyxs
from gdsfactory.technology import LayerLevel, LayerStack, LayerViews, LayerMap
from gdsfactory.generic_tech import get_generic_pdk

from gdsfactory.config import PATH
from gdsfactory.technology import (
    LayerLevel,
    LayerStack,
    LayerView,
    LayerViews,
    LayerMap,
)

Layer = tuple[int, int]

# gf.generic_tech.


class MyLayerMap(LayerMap):
    WG: Layer = (1, 0)
    WGCLAD: Layer = (111, 0)
    BOX: Layer = (112, 0)
    PORT: Layer = (1, 10)
    DESIGN: Layer = (200, 0)
    GUESS: Layer = (201, 0)


LAYER = MyLayerMap

nm = 1e-3
thickness_wg = 220 * nm
thickness_slab_deep_etch = 90 * nm
thickness_slab_shallow_etch = 150 * nm
box_thickness = 0.5
thickness_clad = 0.5

layer_core = LogicalLayer(layer=LAYER.WG)
layer_clad = LogicalLayer(layer=LAYER.WGCLAD)
layer_box = LogicalLayer(layer=LAYER.BOX)
design_layer = LogicalLayer(layer=DESIGN_LAYER)

LAYER_STACK.layers.update(dict(
    core=LayerLevel(
        layer=layer_core,
        thickness=thickness_wg,
        zmin=0,
        material="Si",
        mesh_order=2,
    ),
    box=LayerLevel(
        layer=layer_box,
        thickness=box_thickness,
        zmin=-box_thickness,
        material="SiO2",
        mesh_order=9,
    ),
    clad=LayerLevel(
        layer=layer_clad,
        zmin=0.0,
        material="SiO2",
        thickness=thickness_clad,
        mesh_order=10,
    ),
    design=LayerLevel(
        layer=design_layer,
        zmin=0.0,
        material=None,
        thickness=.0001,
        mesh_order=100,
    ),
))

pdk = get_generic_pdk()
pdk.activate()
LAYER_VIEWS = pdk.layer_views
LAYER_VIEWS.layer_views["WGCLAD"].visible = True
# for l in LAYER_STACK.layers.values():
#     print(l.layer)
#     print(l.layer.layer)
