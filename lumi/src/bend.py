import luminescent as lumi
import gdsfactory as gf
import numpy as np

N = 3  # 3D or 2D
keys = ["2,1"]  # same as keys=["o2@0,o1@0"]
nres = 40

radius = 1.5
wavelengths = 1.55
name = f"bend_R{radius}"

# radius = 2
# wavelengths = lumi.wavelength_range(center=1.35, bandwidth=.4, length=3)
# name = f"bend_R{radius}_multi"

# N = 3  # 2D or 3D
# radius = 5
# wavelengths = 1.55
# name = f"bend_R{radius}"

# N = 3
# radius = 5
# wavelengths = np.linspace(1.5, 1.6, 5)  # number or list or array
# name = f"bend_R{radius}_multi"

N = 2
radius = 5
wavelengths = np.linspace(1.5, 1.6, 5)  # number or list or array
nres = 20
name = f"bend_R{radius}_multi_{N}D"

c = gf.components.bend_circular(radius=radius, allow_min_radius_violation=True)
lumi.write_sparams(c, name=name, wavelengths=wavelengths,
                   nres=nres, keys=keys, N=N, run=False)
# sol = lumi.load_solution(name=name)
