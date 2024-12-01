import luminescent as lumi
import gdsfactory as gf

N = 2

radius = 2
wavelengths = 1.1
name = f"bend_R{radius}__{wavelengths}"
# radius = 2
# wavelengths = [1.1, 1.55]
# name = f"bend_R{radius}_multi"

c = gf.components.bend_circular(radius=radius, allow_min_radius_violation=True)
c.plot()
lumi.write_sparams(c, name=name, wavelengths=wavelengths, keys=["2,1"],  # same as keys=["o2@0,o1@0"]
                   dx=0.05, N=N)
sol = lumi.load_solution(name=name)
