from tidy3d import material_library
from tidy3d.constants import C_0
import tidy3d as td
C0 = td.constants.C_0
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
