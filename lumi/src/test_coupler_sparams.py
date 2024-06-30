import gdsfactory as gf
import json
import pprint as pp
import luminescent as lumi

# import gplugins.luminescent as gl
# from gplugins.luminescent.generic_tech import LAYER, LAYER_STACK
# from gplugins.luminescent.utils import add_bbox
# from gplugins.luminescent.constants import XMARGIN

# c = gf.components.straight(.51)
# # c = gf.components.coupler(gap=0.1, length=4, dx=8, dy=4)
# c = add_bbox(c, layers=[LAYER.WGCLAD, LAYER.BOX],
#              nonport_margin=XMARGIN)
# c.show()

# sol = write_sparams(c, wavelengths=[1.55],
#                     # keys=["1,1", "2,1", "3,1", "4,1"],
#                     dx=0.05,                       approx_2D=False,
#                     gpu=None,                       plot=True)
sol = load_solution()
pp.pprint(sol)
