from pprint import pprint
import luminescent as lumi

name = "1x4_splitter"
c = lumi.gcells.mimo(west=1, east=4, l=6.0, w=6.0, wwg=.5, taper=.05, )
targets = {
    "tparams": {1.55: {"2,1": 0.25, "3,1": 0.25}},
    "phasediff": {1.55: {"2,3": 0.0}},
}

prob = lumi.gcell_problem(
    c, targets, name=name,
    symmetries=[1], lvoid=0.1, lsolid=0.1, dx=0.1,
    N=True, stoploss=.03, iters=40)
sol = lumi.solve(prob)
