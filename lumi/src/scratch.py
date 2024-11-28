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
# print(m, n)
x = np.linspace(0.5*dx, (m-.5)*dx, m)
y = np.linspace(0.5*dx, (n-.5)*dx, n)

i = round(n/2)


def ϵfunc(x_, y_):
    m, n = len(x_), len(y_)
    return resize(eps, (m, n))


tol = 1e-6
# plot = True
# neigs = 2
solver = EMpy.modesolvers.FD.VFDModeSolver(
    λ, x, y, ϵfunc,  "0000").solve(neigs, tol)
# solvers = EMpy.modesolvers.FD.VFDModeSolver(
#     λ, x, y1, ϵfunc1, "SS00").solve(neigs, tol)

# fig = pylab.figure()
# fig.add_subplot(2, 3, 1)
# Hy = np.transpose(solver.modes[-1].get_field("Hy", x, y))
# pylab.imshow(abs(Hy))
# fig.add_subplot(2, 3, 2)
# Hy = np.transpose(solver1D.modes[-1].get_field("Hy", x, y1))
# pylab.imshow(abs(Hy))
# fig.add_subplot(2, 3, 4)
# e = np.transpose(ϵfunc(x, y))
# pylab.imshow(e)
# fig.add_subplot(2, 3, 5)
# e = np.transpose(ϵfunc1(x, y1))
# pylab.imshow(e)
# pylab.show()

modes = [{k: m.get_field(k, x, y) for k in [
    "Ex", "Ey", "Ez", "Hx", "Hy", "Hz"]} for m in solver.modes]
neffs = [np.real(m.neff) for m in solver.modes]
np.savez(os.path.join(path, 'modes.npz'), **modes)
