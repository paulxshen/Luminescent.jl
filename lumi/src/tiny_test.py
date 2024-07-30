from pprint import pprint
import luminescent as lumi

name = "demux"
c = lumi.gcells.mimo(west=1, east=1, l=.7, w=.7,  wwg=.5)
targets = {
    1.55: {
        "2,1": 1.0
    }}
# c.show()

prob = lumi.inverse_design_problem(
    c, tparam_targets=targets,
    lmin=0.2, dx=0.1, maxiters=2, eta=10., approx_2D=True, )
# path="precompile_execution")
# lmin=0.2, dx=0.1, maxiters=2, eta=10., approx_2D=True, gpu="CUDA")
sol = lumi.solve(prob)

lumi.show_solution()
print("post optim tparams:")
pprint(sol["tparams"])
