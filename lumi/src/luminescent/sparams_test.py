from regex import F
import gplugins.luminescent as gl
import gdsfactory as gf

from gplugins.luminescent.generic_tech import LAYER_MAP, LAYER_STACK
from gplugins.luminescent.utils import add_bbox
from gplugins.luminescent.constants import *

c = gf.components.straight(.51)
# c = gf.components.coupler(gap=0.05, length=4, dx=8, dy=4)
c = add_bbox(c, layers=[LAYER_MAP.WGCLAD, LAYER_MAP.BOX],
             nonport_margin=YMARGIN)
c.show()
sol = gl.write_sparams(c,
                       wavelengths=[1.3, 1.55],
                       #    wavelengths=[1.55],
                       keys=["1,2"],
                       layer_stack=LAYER_STACK, dx=0.05,
                       approx_2D=True,
                       #    approx_2D=False,
                       use_gpu=True,
                       plot=True
                       )
