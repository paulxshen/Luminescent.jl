import itertools
from pprint import pprint
from time import sleep
import luminescent as lumi

name = "ubend"
wwg = .5
d = .5
c = lumi.gcells.mimo(west=[wwg/2, 1.5*wwg+d], l=3, w=2*wwg+d,  wwg=wwg)
targets = {
    1.55: {
        "1,2": 1.0
    }}
# c.show()

prob = lumi.inverse_design_problem(
    c, tparam_targets=targets,
    lmin=0.2, dx=0.05, iters=60, eta=4., approx_2D=True, gpu="CUDA", dev=True)
sol = lumi.solve(prob, )
