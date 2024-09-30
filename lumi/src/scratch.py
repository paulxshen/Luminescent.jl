import itertools
from pprint import pprint
import luminescent as lumi
# lumi.make_training_movie(name="mode_converter")
# lumi.make_simulation_movie(name="mode_converter")
c = lumi.gcells.mimo(west=1, east=1, l=1, w=1,  wwg=.5)
targets = {"tparams": {1.55: {"2,1": 1.0}}}
for (approx_2D, gpu, dtype, save_memory) in itertools.product(
    [True,],
    # [None, "CUDA"],
    [None, ],
    ["f32"],
    # ["f32", "f16"],
    [False],
):
    prob = lumi.gcell_problem(
        c, targets,
        # bbox_layer=LAYER.WAFER,
        lvoid=0.2, lsolid=.2, dx=0.1, iters=2,
        approx_2D=approx_2D, gpu=gpu, dtype=dtype, save_memory=save_memory,
        run=False)
    sol = lumi.solve(prob, run=False)
    sleep(1)
