import os
import luminescent as lumi

path = os.path.join("runs", "splitter")
c = lumi.mimo(west=1, east=2, l=4.0, w=2.0, wwg=.5)
targets = {
    "tparams": {1.55: {"2,1": 0.5}},
}

lumi.make_pic_inv_problem(
    path, c, targets,
    nres=20,  symmetries=[1],
    lvoid=0.1, lsolid=.1,
    iters=50, stoploss=.05,
    approx_2D_mode="TE", dtype="float32")  # ),gpu="CUDA")
