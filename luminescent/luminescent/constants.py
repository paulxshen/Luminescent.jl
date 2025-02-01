from gdsfactory.generic_tech import get_generic_pdk
import os
ZMARGIN = .3
XYMARGIN = XMARGIN = YMARGIN = .6

pdk = get_generic_pdk()
pdk.activate()
LAYER_VIEWS = pdk.layer_views
LAYER_VIEWS.layer_views["WGCLAD"].visible = True
DESIGN_LAYER = (1000, 1)

eps0 = 8.854187817e-12
