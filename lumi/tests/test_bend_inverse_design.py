from copy import deepcopy
from pprint import pprint
import gdsfactory as gf
import luminescent as lumi

c = lumi.gcells.bend(r=.5, wwg=.5)
targets = {
    "tparams": {
        1.55: {
            "2,1": 1
        }}}
c.show()

prob = lumi.inverse_design_problem(
    c, targets=targets, wavelengths=[1.55],
    lmin=0.1, dx=0.05, maxiters=20, approx_2D=True)
sol = lumi.solve(prob)
pprint(sol)

c0 = deepcopy(c)
c = lumi.apply_design(c,  sol)
c.write_gds("optimal_bend.gds", "")
c.plot()
