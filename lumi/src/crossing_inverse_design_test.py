from pprint import pprint
import luminescent as lumi

name = "crossing"
L = 4.0
c = lumi.gcells.mimo(1, 1, 1, 1, l=L, w=L, wwg=.5)
targets = {
    1.55: {
        "3,1": 1.0
    }}

prob = lumi.inverse_design_problem(
    c, tparam_targets=targets, symmetries=[0, 1, "diag"], lmin=0.2, dx=0.1,
    minloss=.03, maxiters=50, eta=10., approx_2D=True)
sol = lumi.solve(prob)

# sol = lumi.load_solution()
lumi.show_solution()
print("post optim tparams:")
pprint(sol["tparams"])
