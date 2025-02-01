import luminescent as lumi

margin = 5
path = "ant"
lumi.make_sim_prob(
    # path containing geometry/($material or monitors or sources)/*.STL
    path=path,

    # frequencies=[.1, 1, 2, 4, 8],
    frequencies=[5],
    # any unit (GHz in our case) - 1 unit of simulation time taken as 1 period at this frequency
    center_frequency=5,
    center_wavelength=60,  # same unit as in STL (mm in our case)

    sources=[{
        'mode': {'Jy': 1},  # in local coordinate frame!
        'frame': [  # +z pointing out of port
            [1, 0, 0],
            [0,  1, 0],
            [0, 0, 1]
        ]}],

    monitors=[
        {
            'mode': {'Ey': 1, 'Hx': 1},  # in local frame!
            'frame': [
                [1, 0, 0],
                [0,  1, 0],
                [0, 0, 1]
            ]
        },
        {
            'mode': {'Ey': 1, 'Hx': 1},
            'frame': [
                [-1, 0, 0],
                [0,  1, 0],
                [0, 0, -1]
            ]
        },
    ],

    margins=[[margin, margin, 0], [margin, margin, 0]],  # air margin

    materials={
        'gel': {'epsilon': 50},
        'fr4': {'epsilon': 4.3},
        'PEC': {'epsilon': 1000},
    },

    dx=1,

    Ttrans=1,
    Tss=.1
)

# lumi.solve(path)
