# straight waveguide example 
import luminescent as lumi
from luminescent import LAYER
import gdsfactory as gf
import pprint as pp

c = gf.components.straight(1.0)
c = lumi.add_bbox(c, layers=[LAYER.WGCLAD, LAYER.BOX], nonport_margin=.1)
# c.show()

sol = lumi.write_sparams(c, wavelengths=[1.55], keys=["2,1"],
    #  dx=0.025, approx_2D=True, dtype="float16", gpu=None,)
    dx=0.05, approx_2D=False, gpu="CUDA")
# sol = lumi.load_solution()
lumi.show_solution()
pp.pprint(sol)
