import gdsfactory as gf
from luminescent import LAYER
import luminescent as lumi
import os
os.environ['PATH'] += ':/usr/local/LuminescentCompiled/bin'

c = gf.components.straight(.2)
c = lumi.add_bbox(c, layers=[LAYER.WGCLAD, LAYER.BOX], nonport_margin=.1)
sol = lumi.write_sparams(c, wavelengths=[1.55], keys=[
                         "2,1"], dx=0.1, approx_2D=True, gpu="CUDA")
