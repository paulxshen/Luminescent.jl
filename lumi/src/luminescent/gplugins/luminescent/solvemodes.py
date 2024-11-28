import numpy as np
import EMpy
import os
import sys
from skimage.transform import resize

path = sys.argv[1]
data = np.load(os.path.join(path, "args.npz"))
eps = data["eps"]
dx = data["dx"]
λ = data["wl"]
neigs = data["neigs"]

tol = 1e-4
m, n = eps.shape
x = np.linspace(0.5*dx, (m-.5)*dx, m+1)
y = np.linspace(0.5*dx, (n-.5)*dx, n+1)


def ϵfunc(x_, y_):
    return eps
    # m, n = len(x_), len(y_)
    # return resize(eps, (m, n))


tol = 1e-6
solver = EMpy.modesolvers.FD.VFDModeSolver(
    λ, x, y, ϵfunc,  "0000").solve(neigs, tol)

modes = [{k: m.get_field(k, x, y) for k in [
    "Ex", "Ey", "Ez", "Hx", "Hy", "Hz"]} for m in solver.modes]
neffs = [np.real(m.neff) for m in solver.modes]
for i, mode in enumerate(modes):
    np.savez(os.path.join(path, f'mode{i}.npz'), **modes[i])
