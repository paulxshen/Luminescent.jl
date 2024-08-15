import itertools
from pprint import pprint
from time import sleep
import luminescent as lumi
from gdsfactory.generic_tech import LAYER


name = "demux"
c = lumi.gcells.mimo(west=1, east=1, l=.4, w=.8,  wwg=.5)
targets = {
    1.55: {
        "2,1": 1.0
    }}
# c.show()

prob = lumi.inverse_design_problem(
    c, tparam_targets=targets,
    # bbox_layer=LAYER.WAFER,
    # lmin=0.2, dx=0.1, iters=2, eta=10., approx_2D=True, dev=True)  # gpu="CUDA", dev=True)
    lmin=0.2, dx=0.1, iters=2, eta=10., approx_2D=True, gpu="CUDA", dev=True)
sol = lumi.solve(prob, )
sol = lumi.finetune(2)
c = lumi.apply_design(c, sol)

raise ValueError("stop here")

for (approx_2D, gpu, dtype, ) in itertools.product(
    [True,],
    [None, "CUDA"],
    # [None, ],
    ["f32"],
    # ["f32", "f16"],
):
    prob = lumi.inverse_design_problem(
        c, tparam_targets=targets,
        bbox_layer=LAYER.WAFER,
        lmin=0.2, dx=0.1, iters=2, eta=1., approx_2D=approx_2D, dev=True, gpu=gpu, dtype=dtype, run=False)
    sol = lumi.solve(prob, run=False)
    sleep(1)
# path="precompile_execution")
# lmin=0.2, dx=0.1, iters=2, eta=10., approx_2D=True, gpu="CUDA")

# lumi.show_solution()
# print("post optim tparams:")
# pprint(sol["tparams"])
