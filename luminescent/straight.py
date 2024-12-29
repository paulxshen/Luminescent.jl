import os
import luminescent as lumi
from gdsfactory.generic_tech import LAYER, LAYER_STACK
import gdsfactory as gf
import numpy as np

c = gf.components.straight(length=1, width=0.5, layer=LAYER.WG)
# wavelengths = [1.1, 1.3, 1.5]
wavelengths = np.linspace(1.45, 1.65, 5)
wavelengths = 1.55
path = os.path.join("runs", "straight")
lumi.make_pic_sim_prob(path, c, wavelengths=wavelengths, keys=["o2@1,o1@1"],
                       nres=30, N=3)
# lumi.solve(path)
# sol = lumi.lumi_solution()
