from re import T
import gplugins.luminescent as gl
import gdsfactory as gf

from gplugins.luminescent.generic_tech import LAYER, LAYER_STACK
from gplugins.luminescent.utils import add_bbox

c = gf.components.coupler()
c = gf.components.straight(.6)
# cross_section=gf.CrossSection(     sections=[gf.Section(width=.5, layer=LAYER.WG)]))
xmargin = 0.4
zmargin = 0.2
c = add_bbox(c, layers=[LAYER.WGCLAD, LAYER.BOX],
             nonport_margin=xmargin)
c.show()
sol = gl.write_sparams(
    c, wavelengths=[1.55],
    xmargin=xmargin,
    zmargin=zmargin,
    layer_stack=LAYER_STACK, dx=0.05, approx_2D=True)
