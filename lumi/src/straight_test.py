
import itertools
from time import sleep
import luminescent as lumi
from gdsfactory.generic_tech import LAYER
import gdsfactory as gf
import pprint as pp

c = gf.components.straight(.5)

# sol = lumi.write_sparams(
#     c, wavelengths=[1.55], keys=["2,1"],
#     bbox_layer=LAYER.WAFER,
#     # bbox_layer=[LAYER.WAFER, LAYER.SLAB90],
#     # c, wavelengths=[1.55], keys=["o2@1,o1@1"],
#     #  dx=0.025, approx_2D=True, gpu=None,)
#     dx=0.1, approx_2D=False, gpu="CUDA",)  # dtype="16", dev=True,)
# lumi.show_solution()
# pp.pprint(sol)
# raise ValueError("stop here")

sleep(1)
for (approx_2D, gpu, dtype, wavelengths) in itertools.product(
    [True, False],
    [None, "CUDA"],
    #  ["f32"]
    ["f32", "f16"],
        [[1.55], ]):
    lumi.write_sparams(c, wavelengths=wavelengths, keys=["2,1"], dx=0.1,
                       bbox_layer=LAYER.WAFER,
                       approx_2D=approx_2D, gpu=gpu, dtype=dtype, run=False)
    sleep(1)
#  dx=0.025, approx_2D=True, gpu=None,)
#  path="precompile_execution")  # dtype="float16")
#  dx=0.1, approx_2D=True, gpu="CUDA", dev=True)  # dtype="float16")
# sol = lumi.load_solution()
