import itertools
from pprint import pprint
from time import sleep
import luminescent as lumi
from gdsfactory.generic_tech import LAYER


name = "demux"
c = lumi.gcells.mimo(west=1, east=1, l=1, w=1,  wwg=.5)
targets = {"tparams": {
    1.55: {
        "2,1": 1.0
    }}}
# c.show()

# prob = lumi.gcell_problem(
#     c, targets,
#     # bbox_layer=LAYER.WAFER,
#     # lmin=0.2, dx=0.1, iters=2, eta=10., approx_2D=False, save_memory=True)  # gpu="CUDA", dev=True)
#     lmin=0.2, dx=0.1, iters=2, approx_2D=True, save_memory=False)  # gpu="CUDA", dev=True)
# # lmin=0.2, dx=0.1, iters=2, eta=10., approx_2D=True, gpu="CUDA", dev=True)
# sol = lumi.solve(prob, )
# # sol = lumi.finetune(2)
# sol = lumi.load_solution()
# c = lumi.apply_design(c, sol)
# c.show()
# # raise ValueError("stop here")

for (approx_2D, gpu, dtype, save_memory) in itertools.product(
    [True,],
    # [None, "CUDA"],
    [None, ],
    ["f32"],
    # ["f32", "f16"],
    [True, False],
):
    prob = lumi.gcell_problem(
        c, targets,
        bbox_layer=LAYER.WAFER,
        lmin=0.2, dx=0.1, iters=3,
        approx_2D=approx_2D, gpu=gpu, dtype=dtype, save_memory=save_memory,
        run=False)
    sol = lumi.solve(prob, run=False)
    sleep(1)
