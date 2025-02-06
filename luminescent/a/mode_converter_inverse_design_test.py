# recommended RAM: >16G
import luminescent as lumi

path = "mode_converter"  # can be any string
c = lumi.gcells.mimo(west=1, east=1, l=5.0, w=2.4,
                     wwg=.5, taper=.05)
targets = {"tparams": {1.55: {"o2@1,o1@0": 1.0}}}

prob = lumi.make_pic_inv_problem(
    c, targets, path,
    N=2,  nres=15,
    lvoid=0.15, lsolid=.15,
    iters=100, stoploss=.05, )
# lumi.solve(prob)
lumi.solve(prob, run=False)
