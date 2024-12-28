import itertools
from pprint import pprint
from time import sleep
import luminescent as lumi
from gdsfactory.generic_tech import LAYER


name = "tiny"
# lumi.load_res()
# raise NotImplementedError("This is a stub")
c = lumi.gcells.mimo(west=1, east=1, l=1, w=1,  wwg=.5)
targets = {"tparams": {
    1.55: {
        "2,1": 1.0
    }}}
prob = lumi.make_pic_inv_prob(
    c, targets, name=name,
    lvoid=0.2, dx=0.1, iters=2,
    N=2, gpu=None,)
sol = lumi.solve(prob, run=False)
