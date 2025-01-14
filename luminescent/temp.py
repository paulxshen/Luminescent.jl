# import luminescent as lumi
# import gdsfactory as gf
# import numpy as np
# import os

# radius = 5
# c = gf.components.bend_circular(radius=radius)
# # c.plot()

# path = os.path.join("runs", f"bend_R{radius}")
# wavelengths = 1.55
# keys = ["2,1"]  # same as keys=["o2@0,o1@0"]
# nres = 30
# dtype = "float32"

# lumi.make_pic_sim_prob(path, c, wavelengths=wavelengths,
#                        nres=nres, keys=keys, dtype=dtype)
# sol = lumi.load_res(path)

# import os
# import luminescent as lumi

# path = os.path.join("runs", "splitter")
# c = lumi.mimo(west=1, east=2, l=4.0, w=2.0, wwg=.5, taper=.05, )
# targets = {
#     "tparams": {1.55: {"2,1": 0.5}}  # , "3,1": 0.5}},
# }

# lumi.make_pic_inv_prob(
#     path, c, targets,
#     nres=15,  symmetries=[1],
#     lvoid=0.2, lsolid=.1,
#     iters=50, stoploss=.05,
#     approx_2D_mode="TE")

# import os
# import luminescent as lumi

# path = os.path.join("runs", "demux")
# l = 6.0
# w = 4.0
# wwg = .5
# c = lumi.mimo(west=1, east=[2*wwg, w-2*wwg], l=l, w=w, wwg=wwg, taper=0.05)
# targets = {"tparams": {
#     1.2: {"2,1": 1.0},
#     # 1.8: {"2,1": 1.0},
#     1.8: {"2,1": 1.0, "3,1": .0},
# }}

# lumi.make_pic_inv_prob(
#     path, c, targets,
#     lvoid=0.2, lsolid=0.1, nres=15,
#     approx_2D_mode="TE", Ttrans="1.4x",
#     stoploss=.1, iters=200, dtype="float32")
# sol = lumi.load_res(path)


# 30G RAM
import os
import luminescent as lumi

path = os.path.join("runs", "demux")
c = lumi.mimo(west=1, east=1, south=1, l=4.0, w=4.0, wwg=.5, taper=0.05)
targets = {"tparams": {
    1.2: {"3,1": 1.0},
    1.80: {"2,1": 1.0},
}}

lumi.make_pic_inv_prob(
    path, c, targets,
    lvoid=0.15, lsolid=0.1, nres=30,
    approx_2D_mode="TE", stoploss=.1, iters=200, dtype="float16")  # ,gpu="CUDA")
lumi.solve(path)
