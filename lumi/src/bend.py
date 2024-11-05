import luminescent as lumi
from luminescent import MATERIALS
from gdsfactory.generic_tech import LAYER, LAYER_STACK
import gdsfactory as gf
import pprint as pp

c = gf.components.bend_euler(5)
# c = gf.components.bend_circular(5, allow_min_radius_violation=True)
c.plot()
c.show()
# raise ValueError("This is a test error")
sol = lumi.write_sparams(c, name="bend_euler", wavelengths=1.55, keys=["2,1"],  # same as keys=["o2@0,o1@0"]
                         dx=0.05, N=3, run=False)
# lumi.show_solution()
