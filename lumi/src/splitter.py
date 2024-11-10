from pprint import pprint
import luminescent as lumi

name = "1x2_splitter"
c = lumi.gcells.mimo(west=1, east=2, l=4, w=2.4, wwg=.5, taper=.025, )
targets = {
    "tparams": {1.55: {"2,1": 0.5}},
}

prob = lumi.gcell_problem(
    c, targets, name=name,
    N=2,   symmetries=[1], lvoid=0.2,  dx=0.1,
    eta=.5, iters=3, stoploss=.01, )
sol = lumi.solve(prob)
# name = "1x4_splitter"
# c = lumi.gcells.mimo(west=1, east=4, l=6.0, w=6.0, wwg=.5, taper=.05, )
# targets = {
#     "tparams": {1.55: {"2,1": 0.25, "3,1": 0.25}},
#     "phasediff": {1.55: {"2,3": 0.0}},
# }

# prob = lumi.gcell_problem(
#     c, targets, name=name,
#     symmetries=[1], lvoid=0.1, lsolid=0.1, dx=0.1,
#    N=2, stoploss=.03, iters=40)
# sol = lumi.solve(prob)
