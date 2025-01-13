import numpy as np

try:
    np.finfo(np.float16)
    print("float16 is supported.")
except ValueError:
    print("float16 is not supported.")
