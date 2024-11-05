import matplotlib.pyplot as plt
import numpy as np

# Create some 3D data
x, y, z = np.indices((8, 8, 8))
voxelarray = (x < 3) & (y < 3) & (z < 3)

# Create a figure and axes
fig = plt.figure()
ax = fig.add_subplot(projection='3d')

# Plot the voxels
ax.voxels(voxelarray, edgecolor='k')

plt.show()
