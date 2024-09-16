from pprint import pprint
import luminescent as lumi

# name = "1x2_splitter"
# c = lumi.gcells.mimo(west=1, east=2, l=3.0, w=3.0, wwg=.5, name=name)
# targets = {"tparams": {1.55: {"2,1": 0.5}}}

# prob = lumi.gcell_problem(
#     c, targets,
#     symmetries=[1], lmin=0.15, dx=0.05,
#     approx_2D=True, iters=4,
#     wd="runs")
# sol = lumi.solve(prob)

sol = lumi.load_solution(wd="runs")
sol["optimized_component"].show()
    