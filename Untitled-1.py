import jax
import jax.numpy as jnp


import matplotlib.pyplot as plt
import numpy as onp
# from skimage import measure
# import gifcm

import invrs_gym
# import invrs_opt
challenge = invrs_gym.challenges.ceviche_lightweight_waveguide_bend()
params = challenge.component.init(jax.random.PRNGKey(0))
density = challenge.component.ceviche_model.density(params.array)
1
onp.save("bend.npy",density)