import luminescent as lumi
from luminescent import MATERIALS
from gdsfactory.generic_tech import LAYER, LAYER_STACK
import gdsfactory as gf
import pprint as pp

wg = gf.components.straight(length=2, width=0.5, layer=LAYER.WG)
wg.plot()

c = gf.Component()
c << wg
c.add_ports(wg.ports)
# name="wg_TE0"
sol = lumi.write_sparams(c, wavelength=1.55, keys=["o2@0,o1@0"],  # same as keys=["o2@0,o1@0"]
                         core_layer=LAYER.WG,   bbox_layer=LAYER.WAFER,  # defaults
                         layer_stack=LAYER_STACK, materials=MATERIALS,  # defaults
                         dx=0.1, approx_2D=False, dtype="float32",)  # gpu="CUDA",)  # or gpu=None
lumi.show_solution()
