import numpy as np
import EMpy
import os
import sys

path = sys.argv[1]
data = np.load(os.path.join(path, "args.npz"))
eps = data["eps"]
dl = data["dl"]
λ = data["center_wavelength"]
neigs = data["neigs"]
name = data["name"]

m, n = eps.shape
# print(m, n)
m += 1
n += 1
x = np.linspace(0.5*dl, (m-.5)*dl, m)
y = np.linspace(0.5*dl, (n-.5)*dl, n)


def ϵfunc(x_, y_):
    return eps


tol = 1e-6
solver = EMpy.modesolvers.FD.VFDModeSolver(
    λ, x, y, ϵfunc,  "0000").solve(neigs, tol)

modes = [{k: m.get_field(k, x, y) for k in [
    "Ex", "Ey", "Ez", "Hx", "Hy", "Hz"]} for m in solver.modes]
neffs = [np.real(m.neff) for m in solver.modes]
for i, mode in enumerate(modes):
    # print(mode["Ex"].shape)
    np.savez(os.path.join(path, f'{name}_mode_{i}.npz'), **modes[i])
