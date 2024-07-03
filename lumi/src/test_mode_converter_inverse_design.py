from pprint import pprint
import gdsfactory as gf
import luminescent as lumi

name = "mode_converter"
c = lumi.gcells.mimo(l=3.0, w=2.0, nwest=1, neast=1, wwg=.5)
targets = {
    1.55: {
        "o2@1,o1@0": 1.0
    }}

prob = lumi.inverse_design_problem(
    c, tparam_targets=targets, symmetries=[1],
    lmin=0.2, dx=0.05, maxiters=30, eta=10., approx_2D=True)
sol = lumi.solve(prob)

# sol = lumi.load_solution()
print("post optim tparams:")
pprint(sol["after"]["tparams"])

c = sol["after"]["component"]
# c.write_gds(f"optimal_{name}.gds", "")
