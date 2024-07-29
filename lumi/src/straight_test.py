
import luminescent as lumi
from luminescent import LAYER
import gdsfactory as gf
import pprint as pp

c = gf.components.straight(.5)
c = lumi.add_bbox(c, layers=[LAYER.WGCLAD, LAYER.BOX], nonport_margin=.1)
# c.show()

sol = lumi.write_sparams(c, wavelengths=[1.55], keys=["2,1"],
                         #  dx=0.025, approx_2D=True, gpu=None,)
                         dx=0.1, approx_2D=False, gpu="CUDA", dev=True)  # dtype="float16")
#  dx=0.1, approx_2D=True, gpu="CUDA", dev=True)  # dtype="float16")
# sol = lumi.load_solution()
lumi.show_solution()
pp.pprint(sol)
