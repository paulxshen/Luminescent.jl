# recommended RAM: >16G
from pprint import pprint
import luminescent as lumi

name = "mode_converter"  # can be any string
c = lumi.gcells.mimo(west=1, east=1, l=6.0, w=3.0,
                     wwg=.5, taper=.05)
targets = {"tparams": {1.55: {"o2@1,o1@0": 1.0}}}

prob = lumi.gcell_problem(
    c, targets, name=name,
    lvoid=0.2, dx=0.1,
    N=2, iters=60, stoploss=.03)
lumi.solve(prob)
sol = lumi.load_solution(name)
