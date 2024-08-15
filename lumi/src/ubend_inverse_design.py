from pprint import pprint
import luminescent as lumi

name = "ubend"
wwg = .5
gap = .5
c = lumi.gcells.mimo(west=[wwg/2, 1.5*wwg+gap], l=4.0, w=2*wwg+gap,  wwg=wwg)
targets = {1.55: {"2,1": 1.0}}
# c.show()

prob = lumi.inverse_design_problem(
    c, tparam_targets=targets,
    lmin=0.2, dx=0.05, iters=40, eta=2., approx_2D=True)
sol = lumi.solve(prob)

# sol = lumi.load_solution()
lumi.show_solution()
print("post optim tparams:")
pprint(sol["tparams"])

c = lumi.apply_design(c, sol)
# c.write_gds(f"optimal_{name}.gds", "")
