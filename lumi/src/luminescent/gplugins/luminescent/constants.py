from gdsfactory.generic_tech import get_generic_pdk
ZMARGIN = .3
XYMARGIN = XMARGIN = YMARGIN = .6
PATH = "lumi_runs"

pdk = get_generic_pdk()
pdk.activate()
LAYER_VIEWS = pdk.layer_views
LAYER_VIEWS.layer_views["WGCLAD"].visible = True
DESIGN_LAYER = (1000, 1)
