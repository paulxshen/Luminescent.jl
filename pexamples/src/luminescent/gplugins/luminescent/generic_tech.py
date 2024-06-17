from gdsfactory.generic_tech import LAYER_STACK, LayerMap, get_generic_pdk, LAYER
from gdsfactory.generic_tech import LAYER_STACK, LayerMap
from gdsfactory.technology import (
    LayerLevel,
    LayerStack,
    LayerView,
    LayerViews,
    LayerMap,
)
from gdsfactory.typings import Layer


class MyLayerMap(LayerMap):
    WG: Layer = (1, 0)
    WGCLAD: Layer = (111, 0)
    BOX: Layer = (112, 0)
    PORT: Layer = (1, 10)
    DESIGN: Layer = (200, 0)
    GUESS: Layer = (201, 0)


LAYER_MAP = MyLayerMap()

box_thickness = 0.5
thickness_clad = 0.5
LAYER_STACK.layers.update(dict(
    box=LayerLevel(
        layer=LAYER_MAP.BOX,
        thickness=box_thickness,
        zmin=-box_thickness,
        material="sio2",
        mesh_order=9,
    ),
    clad=LayerLevel(
        layer=LAYER_MAP.WGCLAD,
        zmin=0.0,
        material="sio2",
        thickness=thickness_clad,
        mesh_order=10,
    ),
    design=LayerLevel(
        layer=LAYER_MAP.DESIGN,
        zmin=0.0,
        # material="sio2",
        thickness=.0001,
        mesh_order=100,
    ),
))

pdk = get_generic_pdk()
pdk.activate()
LAYER_VIEWS = pdk.layer_views
LAYER_VIEWS.layer_views["WGCLAD"].visible = True
