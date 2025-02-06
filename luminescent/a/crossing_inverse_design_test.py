from pprint import pprint
import luminescent as lumi

path = "crossing"
L = 4.0
c = lumi.gcells.mimo(1, 1, 1, 1, l=L, w=L, wwg=.5)
targets = {
    1.55: {
        "3,1": 1.0
    }}

prob = lumi.make_pic_inv_problem(
    c, tparam_targets=targets, symmetries=[0, 1, "diag"], lmin=0.2, dx=0.1,
    stoploss=.03, iters=50, eta=10., N=2)
sol = lumi.solve(prob)

# sol = lumi.load_solution()
lumi.show_solution()
print("post optim tparams:")
pprint(sol["tparams"])
