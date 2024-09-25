from pprint import pprint
import luminescent as lumi

name = "1x2_splitter"
c = lumi.gcells.mimo(west=1, east=2, l=4.0, w=2.0, wwg=.5, name=name)
targets = {"tparams": {1.55: {"2,1": 0.5, "3,1": 0.5}}}

prob = lumi.gcell_problem(
    c, targets,
    symmetries=[1], lvoid=0.1, dx=0.05,
    approx_2D=True, iters=30, stoploss=.03)
sol = lumi.solve(prob)
