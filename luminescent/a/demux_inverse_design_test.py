# RAM: 32G
from pprint import pprint
import luminescent as lumi

path = "demux"
c = lumi.gcells.mimo(west=1, east=2, l=5.0, w=3.0, wwg=.5, path)
targets = {"tparams": {
    1.55: {"2,1": 1.0},
    1.10: {"3,1": 1.0},
}}

prob = lumi.make_pic_inv_problem(
    c, targets,
    lvoid=0.1, lsolid=0.1, dx=0.1,
    N=2, iters=50)
sol = lumi.solve(prob)
