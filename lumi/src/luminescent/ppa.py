# from gcells import ubend
from operator import inv
from cv2 import solve
from regex import E
from gcells import *
from prob import *
import pcells
import gdsfactory as gf
# import picToGDS as pg
import luminescent.gplugins.luminescent.utils as utils
from generic_tech import LAYER_STACK, LAYER, LAYER_VIEWS

dx = .05
wwg = 0.4
r = .4
h = 0.2
l = 3
lmin = 0.1


c, sources, monitors, design = ubend_template(
    r,  l, wwg, dx=dx, wavelengths=1.55, LAYER=LAYER)
c = add_bbox(c, layers=[LAYER.BOX, LAYER.WGCLAD], margin=.2)
prob = gcell_problem(
    c, lmin=lmin, dx=dx, sources=sources, monitors=monitors, design=design, layer_stack=LAYER_STACK)
solve(prob)
raise Exception("Not implemented")
utils.pic2gds("opt/ubend.png", dx)
init = gf.import_gds("image.gds", "TOP",  read_metadata=True).rotate(90)
ubend = ubend_template(r, wwg, h, l, dx=dx, dir="opt", lmin=lmin, init=init)
ubend.show()
ps = pcells.phase_shifter(l=20, wwg=wwg, nbends=4, ubend=ubend)
res = pcells.ubend_resonator(l=50, wwg=wwg, ubend=ubend)
ps.show()
res.show()
