# recommended RAM: >16G
from pprint import pprint
import luminescent as lumi

name = "mode_converter"  # can be any string
c = lumi.gcells.mimo(west=1, east=1, l=4.0, w=2.0, wwg=.5, name=name)
targets = {"tparams": {1.55: {"o2@0,o1@1": 1.0}}}

prob = lumi.gcell_problem(
    c, targets,
    lmin=0.1, dx=0.05,
    approx_2D=True, iters=40, stoploss=.03)
sol = lumi.solve(prob)
