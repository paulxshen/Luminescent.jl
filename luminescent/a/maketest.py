import itertools
import os
import shutil
from pprint import pprint
from time import sleep
import luminescent as lumi
from gdsfactory.generic_tech import LAYER
import gdsfactory as gf

BUILD_RUNS = os.path.join("build", "precompile_execution")

# for filename in os.listdir(BUILD_RUNS):
#     file_path = os.path.join(BUILD_RUNS, filename)
#     if os.path.isfile(file_path) or os.path.islink(file_path):
#         os.unlink(file_path)
#     elif os.path.isdir(file_path):
#         shutil.rmtree(file_path)

nres = 15
c = lumi.gcells.mimo(west=1, east=1, l=1, w=1,  wwg=.4)
targets = {"tparams": {1.55: {"2,1": 1.0}}}
for (gpu, dtype, ) in itertools.product(
    [None, ],
    ["f32"],
    # ["f32", "f16"],
):
    prob = lumi.make_pic_inv_problem(
        c, targets, name="invtest",
        lvoid=0.2, lsolid=.2,  nres=nres, iters=2,
        N=2, gpu=gpu, dtype=dtype,)  # wd=BUILD_RUNS)
    sol = lumi.solve(prob, run=False)
    sleep(1)
    raise Exception("stop")
c = gf.components.straight(.5,)
i = 1
for (N, gpu, dtype, wavelengths) in itertools.product(
    [2, 3],
    [None, ],
    ["f32"],
        [[1.55], ]):
    lumi.make_pic_sim_problem(c, name=f"simtest{i}",
                              wavelengths=wavelengths, keys=["2,1"], nres=nres,
                              N=N, gpu=gpu, dtype=dtype,
                              wd=BUILD_RUNS, run=False)
    i += 1
    sleep(1)
