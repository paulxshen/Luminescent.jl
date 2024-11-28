import luminescent as lumi
from luminescent import MATERIALS
from gdsfactory.generic_tech import LAYER, LAYER_STACK
import gdsfactory as gf
import pprint as pp

c = gf.components.straight(length=1, width=0.5, layer=LAYER.WG)
lumi.write_sparams(c, wavelengths=1.55, keys=["o2@0,o1@0"],  # same as keys=["o2@0,o1@0"]
                   name="straight",
                   #  core_layer=LAYER.WG,   bbox_layer=LAYER.WAFER,  # defaults
                        #  layer_stack=LAYER_STACK, materials=MATERIALS,  # defaults
                   dx=0.1, N=2)  # ,  gpu="CUDA",)  # or gpu=None
sol = lumi.lumi_solution()
