from pprint import pprint
import gdsfactory as gf
import luminescent as lumi

name = "demux"
c = lumi.gcells.mimo(l=4.0, w=4.0, nwest=1, neast=2, wwg=.5)
targets = {
    1.55: {
        "2,1": 1.0
    },
    .85: {
        "3,1": 1.0
    }}
c.show()

prob = lumi.inverse_design_problem(
    c, tparam_targets=targets,
    lmin=0.2, dx=0.05, maxiters=40, eta=10., approx_2D=True)
sol = lumi.solve(prob)

# sol = lumi.load_solution()
print("post optim tparams:")
pprint(sol["after"]["tparams"])

# c = sol["after"]["component"]
# c.write_gds(f"optimal_{name}.gds", "")
