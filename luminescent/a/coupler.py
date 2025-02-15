import luminescent as lumi
from luminescent import MATERIALS
from gdsfactory.generic_tech import LAYER, LAYER_STACK
import gdsfactory as gf
import pprint as pp

c = gf.components.coupler(length=8, gap=.06, dx=10, dy=5)
c.plot()
sol = lumi.make_pic_sim_problem(c, name="coupler2D", wavelengths=1.55, keys=["2,1", "3,1", "4,1"],  # same as keys=["o2@0,o1@0"]
                                core_layer=LAYER.WG,   bbox_layer=LAYER.WAFER,  # defaults
                                layer_stack=SOI, materials=MATERIALS,  # defaults
                                dx=0.1, N=2, N=3)
# lumi.show_solution()
