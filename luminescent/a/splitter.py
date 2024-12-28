from pprint import pprint
import luminescent as lumi

path = "splitter"
c = lumi.gcells.mimo(west=1, east=2, l=4.0, w=2.0, wwg=.5, taper=.05, )
targets = {
    "tparams": {1.55: {"2,1": 0.5}},
}

prob = lumi.make_pic_inv_prob(
    c, targets, path,
    N=2,  nres=15,  symmetries=[1],
    lvoid=0.15, lsolid=.15,
    iters=50, stoploss=.05, )
lumi.solve(prob, run=False)
# path = "1x4_splitter"
# c = lumi.gcells.mimo(west=1, east=4, l=6.0, w=6.0, wwg=.5, taper=.05, )
# targets = {
#     "tparams": {1.55: {"2,1": 0.25, "3,1": 0.25}},
#     "phasediff": {1.55: {"2,3": 0.0}},
# }

# prob = lumi.make_pic_inv_prob(
#     c, targets, path,
#     symmetries=[1], lvoid=0.1, lsolid=0.1, dx=0.1,
#    N=2, stoploss=.03, iters=40)
# sol = lumi.solve(prob)

# apt install libgl1-mesa-glx
# sudo apt-get install libcairo2-dev
# sudo apt install libxcb-cursor0
