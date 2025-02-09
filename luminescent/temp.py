# 60G RAM
import os
import luminescent as lumi

path = os.path.join("runs", "demux")
c = lumi.mimo(west=1, east=2, l=5.0, w=3.0, wwg=.5, taper=0.05)
targets = {"tparams": {
    1.2: {"3,1": 1.0},
    1.8: {"2,1": 1.0},
}}

lumi.make_pic_inv_problem(
    path, c, targets,
    lvoid=0.15, lsolid=0.1, nres=30,
    approx_2D_mode="TE", stoploss=.03, iters=200, dtype="float16")  # ,gpu="CUDA")
lumi.solve(path)
