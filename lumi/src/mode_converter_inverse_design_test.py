# recommended RAM: >16G
from pprint import pprint
import luminescent as lumi

name = "mode_converter"  # can be any string
c = lumi.gcells.mimo(west=1, east=1, l=6.0, w=3.0,
                     wwg=.5, taper=.05, name=name)
targets = {"tparams": {1.55: {"o2@0,o1@1": 1.0}}}

prob = lumi.gcell_problem(
    c, targets,
    lvoid=0.2, lsolid=.1, dx=0.05,
    approx_2D=True, iters=60, stoploss=.03)
sol = lumi.solve(prob)
