# from gcells import ubend
from gcells import ubend
import pcells
import gdsfactory as gf
# import picToGDS as pg
import utils

dx = .05
# ubend = ubend(.4, .5, .22, 2, dir="opt", )
utils.pic2gds("ubend.png", dx)
init = gf.import_gds("image.gds", "TOP",  read_metadata=True)
ubend = ubend(.4, .5, .22, 2, dir="opt", init=init)
ubend.show()
ps = pcells.phase_shifter(l=10, wwg=0.5, nbends=3, ubend=ubend)
res = pcells.ubend_resonator(l=20, wwg=0.5, ubend=ubend)
ps.show()
res.show()
