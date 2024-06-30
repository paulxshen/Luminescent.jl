from copy import deepcopy
import json
from pprint import pprint
import gdsfactory as gf
import luminescent as lumi
# from luminescent import *
# import luminescent as


# import gplugins.luminescent as gl
# from gplugins.luminescent.generic_tech import LAYER, LAYER_STACK
# from gplugins.luminescent.utils import add_bbox
# from gplugins.luminescent.gcells import bend
# from lumi.src.luminescent.gplugins.luminescent.inverse_design import apply_design

# c = lumi.gcells.bend(r=.8, wwg=.4)
c = lumi.gcells.bend(r=.5, wwg=.5)
targets = {
    "tparams": {
        1.55: {
            "2,1": 1
        }}}
c.show()

prob = lumi.inverse_design_problem(
    c, targets=targets, wavelengths=[1.55],
    lmin=0.1, dx=0.05, maxiters=8, approx_2D=True)
sol = lumi.solve(prob)
print("pre optim tparams:")
print(json.dumps(sol["before"]["tparams"]))
print("post optim tparams:")
print(json.dumps(sol["after"]["tparams"]))

c0 = deepcopy(c)
c = lumi.apply_design(c,  sol)
c.write_gds("optimal_bend.gds", "")
c.plot()
