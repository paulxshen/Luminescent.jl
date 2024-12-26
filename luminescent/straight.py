import luminescent as lumi
from gdsfactory.generic_tech import LAYER, LAYER_STACK
import gdsfactory as gf
import numpy as np

c = gf.components.straight(length=1, width=0.5, layer=LAYER.WG)
# wavelengths = [1.1, 1.3, 1.5]
wavelengths = np.linspace(1.45, 1.65, 5)
wavelengths = 1.55
lumi.write_sparams(c, wavelengths=wavelengths, keys=["o2@0,o1@0"],  # same as keys=["o2@0,o1@0"]
                   name="straight",
                   #  core_layer=LAYER.WG,   bbox_layer=LAYER.WAFER,  # defaults
                        #  layer_stack=LAYER_STACK,   # defaults
                   nres=15, N=2)  # , run=False)
#     dx=0.05, N=3,)  # run=False)  # ,  gpu="CUDA",)  # or gpu=None
# sol = lumi.lumi_solution()
