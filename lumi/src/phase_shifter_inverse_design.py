from pprint import pprint
import luminescent as lumi

name = "phase_shifter"
wwg = .5
c = lumi.gcells.mimo(west=1, east=1, l=5.0, w=5,  wwg=wwg)
# c = lumi.gcells.mimo(west=1, east=1, l=1, w=.5,  wwg=wwg)
# c.show()

prob = lumi.inverse_design_problem(
    c, preset={"name": "phase_shifter", "wavelength": 1.55},
    lmin=0.15, dx=0.05, iters=40, approx_2D=True)
sol = lumi.solve(prob)

# sol = lumi.load_solution()
lumi.show_solution()
print("post optim tparams:")
pprint(sol["tparams"])

c = lumi.apply_design(c, sol)
# c.write_gds(f"optimal_{name}.gds", "")
