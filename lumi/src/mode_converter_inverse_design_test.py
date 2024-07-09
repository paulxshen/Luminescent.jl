from pprint import pprint
import luminescent as lumi

name = "mode_converter"
c = lumi.gcells.mimo(west=1, east=1, l=6.0, w=4.0, wwg=.5)
targets = {
    1.55: {
        "o2@1,o1@0": 1.0
    }}

prob = lumi.inverse_design_problem(
    c, tparam_targets=targets, symmetries=[], lmin=0.2, dx=0.05,
    minloss=.03, maxiters=50, eta=10., approx_2D=True)
sol = lumi.solve(prob)

# sol = lumi.load_solution()
lumi.show_solution(sol)
print("post optim tparams:")
pprint(sol["tparams"])