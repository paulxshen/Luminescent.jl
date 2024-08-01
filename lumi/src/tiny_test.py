from pprint import pprint
from time import sleep
import luminescent as lumi

name = "demux"
c = lumi.gcells.mimo(west=1, east=1, l=.2, w=.7,  wwg=.5)
targets = {
    1.55: {
        "2,1": 1.0
    }}
# c.show()

prob = lumi.inverse_design_problem(
    c, tparam_targets=targets,
    # lmin=0.2, dx=0.1, maxiters=2, eta=10., approx_2D=True, dev=True)  # gpu="CUDA", dev=True)
    lmin=0.2, dx=0.1, maxiters=2, eta=10., approx_2D=True, dev=True, gpu="CUDA", run=False)
sol = lumi.solve(prob, run=False)
sleep(1)
prob = lumi.inverse_design_problem(
    c, tparam_targets=targets,
    lmin=0.2, dx=0.1, maxiters=2, eta=10., approx_2D=True, dev=True,
    run=False)  # gpu="CUDA", dev=True)
# path="precompile_execution")
# lmin=0.2, dx=0.1, maxiters=2, eta=10., approx_2D=True, gpu="CUDA")
sol = lumi.solve(prob, run=False)

lumi.show_solution()
print("post optim tparams:")
pprint(sol["tparams"])
