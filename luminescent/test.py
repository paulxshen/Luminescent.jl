import itertools
import os
import luminescent as lumi
import gdsfactory as gf

# dir = "runs"
dir = os.path.join("build", "precompile_execution")

c = gf.components.straight(length=.1, width=0.5,)
wavelengths = 1.5

path = "a"
approx_2D_mode = "TE"
lumi.make_pic_sim_problem(path, c, wavelengths=wavelengths, keys=[
    "2,1"], nres=15, approx_2D_mode=approx_2D_mode, gpu="CUDA")
lumi.solve(path)
