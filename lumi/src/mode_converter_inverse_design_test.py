# recommended RAM: >16G
from pprint import pprint
import luminescent as lumi

name = "mode_converter"  # can be any string
c = lumi.gcells.mimo(west=1, east=1, l=3.0, w=3.0, wwg=.5, name=name)
targets = {"tparams": {1.55: {"o2@1,o1@0": 1.0}}}

prob = lumi.gcell_problem(
    c, targets,
    lmin=0.15, dx=0.05,
    approx_2D=True, iters=50)
sol = lumi.solve(prob)
