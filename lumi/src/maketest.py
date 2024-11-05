import itertools
import os
import shutil
from pprint import pprint
from time import sleep
import luminescent as lumi
from gdsfactory.generic_tech import LAYER
import gdsfactory as gf

# BUILD_RUNS = os.path.join("build", "precompile_execution")
# for filename in os.listdir(BUILD_RUNS):
#     file_path = os.path.join(BUILD_RUNS, filename)
#     if os.path.isfile(file_path) or os.path.islink(file_path):
#         os.unlink(file_path)
#     elif os.path.isdir(file_path):
#         shutil.rmtree(file_path)


c = lumi.gcells.mimo(west=1, east=1, l=1, w=1,  wwg=.4)
targets = {"tparams": {1.55: {"2,1": 1.0}}}
for (N, gpu, dtype, save_memory) in itertools.product(
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
        N=N, gpu=gpu, dtype=dtype, save_memory=save_memory,
        N=3)
    # N=3, wd=BUILD_RUNS)
    sol = lumi.solve(prob, N=3)
    sleep(1)

# c = gf.components.straight(.5,)
# i = 1
# for (N, gpu, dtype, wavelengths) in itertools.product(
#     [False, True],
#     [None, "CUDA"],
#     ["f32"],
#         [[1.55], ]):
#     lumi.write_sparams(c, name=f"{i}",
#                        wavelengths=wavelengths, keys=["2,1"], dx=0.1,
#                        N=N, gpu=gpu, dtype=dtype,
#                        N=3, wd=BUILD_RUNS)
#     i += 1
#     sleep(1)
