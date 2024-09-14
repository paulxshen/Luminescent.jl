from pprint import pprint
import luminescent as lumi

name = "mode_converter"
c = lumi.gcells.mimo(west=1, east=1, l=5.0, w=4.0, wwg=.5)
targets = {
    1.55: {
        "o2@1,o1@0": 1.0
    }}

prob = lumi.gcell_problem(
    c, tparam_targets=targets, symmetries=[], lmin=0.2, dx=0.05,
    iters=100, eta=4.0, approx_2D=True)
sol = lumi.solve(prob)

# sol = lumi.load_solution()
lumi.show_solution()
print("post optim tparams:")
pprint(sol["tparams"])
