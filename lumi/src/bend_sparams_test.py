import luminescent as lumi
from luminescent import LAYER
import gdsfactory as gf
import pprint as pp

# name="wg_multi"
wg = gf.components.straight(length=0.5, width=0.5, layer=LAYER.WG)

c = gf.Component()
c << wg
c.add_ports(wg.ports)

sol = lumi.write_sparams(c, wavelength=[1.25, 1.55], keys=["2,1"],
                         dx=0.05, approx_2D=False, dtype="float32", gpu="CUDA",)  # or gpu=None
# sol = lumi.load_solution()
lumi.show_solution()
