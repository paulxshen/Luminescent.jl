

class Si():
    epsilon = 3.4757**2


class SiO2():
    epsilon = 1.444**2


class SiN():
    epsilon = 2.0**2


class Ge():
    epsilon = 4.2**2


MATERIALS = {
    "si": Si,
    "sio2": SiO2,
    "Si": Si,
    "SiO2": SiO2,
    "SiN": SiN,
    "sin": SiN,
    "Ge": Ge,
    "ge": Ge,

}
# Ge["epsilon"]
