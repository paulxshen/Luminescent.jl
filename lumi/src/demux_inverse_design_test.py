from pprint import pprint
import luminescent as lumi

name = "demux"
c = lumi.gcells.mimo(west=1, east=2, l=4.0, w=4.0,  wwg=.5)
targets = {
    1.55: {
        "2,1": 1.0
    },
    .85: {
        "3,1": 1.0
    }}
# c.show()

prob = lumi.inverse_design_problem(
    c, tparam_targets=targets,
    lmin=0.2, dx=0.05, maxiters=40, eta=10., approx_2D=True)
sol = lumi.solve(prob)

# sol = lumi.load_solution()
lumi.show_solution(sol)
print("post optim tparams:")
pprint(sol["tparams"])

# c = sol["component"]
# c.write_gds(f"optimal_{name}.gds", "")
