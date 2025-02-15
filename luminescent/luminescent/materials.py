import copy
from gdsfactory.generic_tech import LAYER_STACK, LAYER

MATKEYS = {
    "si": "cSi",
    "Si": "cSi",

    "sio2": "SiO2",
    "sin": "SiN",
    "ge": "Ge",

}


def matname(k):
    if k in MATKEYS:
        return MATKEYS[k]
    return k.capitalize()


MATERIALS = {
    "cSi": {
        'epsilon': 3.48**2},
    "SiO2": {'epsilon': 1.44**2},
    "SiN": {'epsilon': 2.0**2},
    "Ge": {'epsilon': 4.0**2},
    "Si": {'epsilon': 3.48**2},
    'PEC': {'epsilon': 10000},
}

ks = copy.deepcopy(list(MATERIALS.keys()))
for k in ks:
    MATERIALS[k.lower()] = MATERIALS[k]

SOI = {
    'layers': dict()
}
SOI['layers']['core'] = LAYER_STACK['core']
SOI['layers']['default'] = {
    'material': 'SiO2'
}
BBOX_LAYER = (8888, 8888)
