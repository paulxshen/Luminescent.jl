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
