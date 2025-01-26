
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
    'PEC': {'epsilon': 1000},
}
