
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
