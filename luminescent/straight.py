import os
import luminescent as lumi
from gdsfactory.generic_tech import LAYER, LAYER_STACK
import gdsfactory as gf
import numpy as np

c = gf.components.straight(length=2, width=0.5, layer=LAYER.WG)
# wavelengths = [1.1, 1.3, 1.5]
wavelengths = [1.2, 1.8]
path = os.path.join("runs", "straight")
lumi.make_pic_sim_problem(path, c, wavelengths=wavelengths, keys=[
    "2,1"],                       nres=30, approx_2D_mode="TE")
# lumi.solve(path)
# sol = lumi.lumi_solution()
