import luminescent as lumi
from luminescent import eps0
import numpy as np

path = "wg"
center_frequency = 5
center_wavelength = 60
frequencies = [2, 3, 4, 5, 6, 7, 8]
sigma = 1/(center_frequency*1e9)/eps0
Z = 50*(center_frequency*1e9)*eps0*1e3/center_wavelength
print(f"sigma: {sigma}, Z: {Z}")

margin = 10
mode = {'Ey': 1, 'Hx': -1/Z}
# +z pointing out of port
frame2 = [[1, 0, 0],
          [0,  1, 0],
          [0, 0, 1]]
frame1 = [[-1, 0, 0],
          [0,  1, 0],
          [0, 0, -1]]
margins = [[margin, margin, 0], [margin, margin, 0]]  # air margin

materials = {
    'gel': {'epsilon': 50, 'sigma': sigma},
    'fr4': {'epsilon': 4.3, },
    'PEC': {'epsilon': 1000},
}

dx = .8
Ttrans = 3

lumi.make_sim_prob(
    path=path,
    frequencies=frequencies, center_frequency=center_frequency, center_wavelength=center_wavelength,


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
    ],

    materials=materials,
    margins=margins,      dx=dx,    Ttrans=Ttrans,
    # Tss=1
)
