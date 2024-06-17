import gplugins.luminescent as gl
import gdsfactory as gf
# import luminescent as


from gplugins.luminescent.generic_tech import LAYER_MAP, LAYER_STACK
from gplugins.luminescent.utils import add_bbox
from gplugins.luminescent.gcells import bend

# c = gl.gcells.bend(r=.8, wwg=.4)
c = bend(r=.8, wwg=.4)
sparam_targets = {
    1.55: {
        "2,1": {
            "mag": 1.,
            "phase": None
        }}
}
c.show()
prob = gl.inverse_design_problem(
    c, sparam_targets=sparam_targets, wavelengths=[1.55],
    lmin=0.1, dx=0.05, maxiters=25, approx_2D=True)
sol = gl.solve(prob)
