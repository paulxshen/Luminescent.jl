import luminescent as lumi
from luminescent import LAYER
import gdsfactory as gf
import pprint as pp

c = gf.components.bend_euler(radius=.75, allow_min_radius_violation=True)
c = lumi.add_bbox(c, layers=[LAYER.WGCLAD, LAYER.BOX], nonport_margin=.1)
# c.show()

sol = lumi.write_sparams(c, wavelengths=[1.55], keys=["2,1"],
                         dx=0.05, approx_2D=False, gpu=None,)
# sol = lumi.load_solution()
lumi.show_solution()
pp.pprint(sol)
