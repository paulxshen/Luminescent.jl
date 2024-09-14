from pprint import pprint
from time import sleep
import luminescent as lumi
from gdsfactory.generic_tech import LAYER
import gdsfactory as gf

# lumi.show_solution()
# raise Exception("stop here")

# c = gf.components.straight(.5)
# sol = lumi.write_sparams(
#     c, name="straight_waveguide",
#     wavelength=1.55, keys=["2,1"],
#     # c, wavelengths=[1.55], keys=["o2@1,o1@1"],
#     bbox_layer=LAYER.WAFER,
#     # bbox_layer=[LAYER.WAFER, LAYER.SLAB90],
#     dx=0.1, approx_2D=True, gpu="CUDA",)
# # dx=0.1, approx_2D=False, gpu="CUDA",)  # dtype="16", dev=True,)
# lumi.show_solution()

# c = lumi.gcells.mimo(west=1, east=1, l=1, w=1,  wwg=.5)
# targets = {"tparams": {
#     1.55: {
#         "2,1": 1.0
#     }}}
# prob = lumi.inverse_design_problem(
#     c, targets, name="tiny_mimo",
#     lmin=0.2, dx=0.1, iters=2, approx_2D=True, save_memory=False)  # gpu="CUDA", dev=True)
# # lmin=0.2, dx=0.1, iters=2,  approx_2D=True, gpu="CUDA", dev=True)
# sol = lumi.solve(prob, )
# sol = lumi.finetune(2)
lumi.show_solution()
sol = lumi.load_solution()
# pprint(sol)
