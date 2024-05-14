# from gcells import ubend
from regex import E
from gcells import ubend
import pcells
import gdsfactory as gf
# import picToGDS as pg
import utils

dx = .05
wwg = 0.4
# ubend = ubend(.35, .4, .2, 2.5, dx=.025, dir="opt", lmin=0.075)

# raise Exception("Not implemented")
utils.pic2gds("opt/ubend.png", dx)
init = gf.import_gds("image.gds", "TOP",  read_metadata=True).rotate(90)
ubend = ubend(.4, wwg, .2, 3, dx=dx, dir="opt", lmin=0.1, init=init)
ubend.show()
ps = pcells.phase_shifter(l=20, wwg=wwg, nbends=4, ubend=ubend)
res = pcells.ubend_resonator(l=50, wwg=wwg, ubend=ubend)
ps.show()
res.show()
