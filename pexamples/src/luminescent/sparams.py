from re import T
import gplugins.luminescent as gl
import gdsfactory as gf

from gplugins.luminescent.generic_tech import LAYER_MAP, LAYER_STACK
from gplugins.luminescent.utils import add_bbox

c = gf.components.coupler()
c = gf.components.straight(.6)
# cross_section=gf.CrossSection(     sections=[gf.Section(width=.5, layer=LAYER_MAP.WG)]))
waveguide_width_margin = 0.4
waveguide_height_margin = 0.2
c = add_bbox(c, layers=[LAYER_MAP.WGCLAD, LAYER_MAP.BOX],
             nonport_margin=waveguide_width_margin)
c.show()
sol = gl.write_sparams(
    c, wavelengths=[1.55],
    waveguide_width_margin=waveguide_width_margin,
    waveguide_height_margin=waveguide_height_margin,
    layer_stack=LAYER_STACK, dx=0.05, approx_2D=True)
