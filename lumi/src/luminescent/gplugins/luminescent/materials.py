

class Si():
    epsilon = 3.4757**2


class SiO2():
    epsilon = 1.444**2


MATERIAL_LIBRARY = {
    "si": Si,
    "sio2": SiO2,
    "Si": Si,
    "SiO2": SiO2,
}
