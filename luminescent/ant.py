import luminescent as lumi
from luminescent import eps0
import numpy as np

path = "drive/MyDrive/reflex/ant"
center_frequency = 5
center_wavelength = 60
frequencies = [2, 3, 4, 5, 6, 7, 8]
sigma = 1/(center_frequency*1e9)/eps0
Z = 50*(center_frequency*1e9)*eps0*1e3/center_wavelength
print(f"sigma: {sigma}, Z: {Z}")

margin = 30
mode = {'Ey': 1, 'Hx': -1/Z}
# +z pointing out of port
frame1 = [[1, 0, 0],
          [0,  1, 0],
          [0, 0, 1]]
frame2 = [[-1, 0, 0],
          [0,  1, 0],
          [0, 0, -1]]
margins = [[margin, margin, 0], [margin, margin, 0]]  # air margin

materials = {
    'gel': {'epsilon': 50, 'sigma': sigma},
    'fr4': {'epsilon': 4.3, },
    'PEC': {'epsilon': 1000},
}

dx = .8
Ttrans = (3*120*7+60*2)/center_wavelength

lumi.make_sim_prob(
    # path containing geometry/($material or monitors or sources)/*.STL
    path=path,
    frequencies=frequencies, center_frequency=center_frequency, center_wavelength=center_wavelength,
    # frequencies=np.linspace(1, 5, 5).tolist(),
    # any unit (GHz in our case) - 1 unit of simulation time taken as 1 period at this frequency

    # same unit as in STL (mm in our case)


    sources=[{
        'mode': {'Jy': 1},  # in local coordinate frame!
        'frame': frame1}],

    monitors=[
        {
            'mode': mode,  # in local frame!
            'frame': frame1
        },
        {
            'mode': mode,
            'frame': frame2
        },
        {
            'mode': mode,
            'frame': frame1
        },
    ],

    materials=materials,
    margins=margins,      dx=dx,    Ttrans=Ttrans,
    # Tss=1
)

lumi.solve(path)
