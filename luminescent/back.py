import os
import luminescent as lumi


path = os.path.join("runs", "back")
# lumi.load_res()
# raise NotImplementedError("This is a stub")
c = lumi.mimo(west=1, east=1, l=.5, w=.5,  wwg=.5)
targets = {"tparams": {
    1.55: {
        "2,1": 1.0
    }}}
lumi.make_pic_inv_prob(
    path,  c, targets,
    lvoid=0.2, iters=2, nres=15,
    approx_2D_mode="TE", gpu="CUDA",)
lumi.solve(path, dev=True)
