from pprint import pprint
import luminescent as lumi

c = lumi.gcells.mimo(west=1, east=1, l=3.0, w=3.0,
                     wwg=.5, name="mode_converter")
targets = {"tparams": {1.55: {"o2@0,o1@1": 1.0}}}

prob = lumi.gcell_problem(
    c, targets,
    lmin=0.15, dx=0.05,
    approx_2D=True, iters=40)
sol = lumi.solve(prob)
