import luminescent as lumi
import gdsfactory as gf

radius = 5
N = 3
c = gf.components.bend_circular(radius=radius, allow_min_radius_violation=True)
c.plot()

name = f"bend_R{radius}_{N}D"
lumi.write_sparams(c, name=name, wavelengths=1.55, keys=["2,1"],  # same as keys=["o2@0,o1@0"]
                   dx=0.05, N=N)
sol = lumi.load_solution(name=name)
