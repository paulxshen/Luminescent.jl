import os
import luminescent as lumi
from gdsfactory.generic_tech import LAYER, LAYER_STACK
import gdsfactory as gf
import numpy as np

c = gf.components.straight(length=.1, width=0.5, layer=LAYER.WG)
wavelengths = 1.55
path = os.path.join("runs", "tiny")
lumi.make_pic_sim_prob(path, c, wavelengths=wavelengths, keys=[
                       "2,1"],                       nres=15, approx_2D_mode="TE", gpu="CUDA")  # approx_2D_mode="TE")
path = os.path.join("runs", "tiny3")
lumi.make_pic_sim_prob(path, c, wavelengths=wavelengths, keys=[
                       "2,1"],                       nres=15, gpu="CUDA")  # approx_2D_mode="TE")
# lumi.solve(path)
# sol = lumi.lumi_solution()
