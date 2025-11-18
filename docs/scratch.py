import matplotlib.pyplot as plt
import numpy as np

# Create the mesh in polar coordinates
r = np.linspace(0, 1, 100)  # Radial distance
theta = np.linspace(0, 2 * np.pi, 100)  # Angle
R, THETA = np.meshgrid(r, theta)

# Define a function in terms of polar coordinates (e.g., a spiral surface)
Z = np.sin(R * 5) + np.cos(THETA * 3)

# Convert polar coordinates to Cartesian coordinates for plotting
X = R * np.cos(THETA)
Y = R * np.sin(THETA)

# Create the 3D plot
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection="3d")

# Plot the surface
ax.plot_surface(X, Y, Z, cmap="viridis")

# Set labels
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
ax.set_title("3D Surface Plot in Polar Coordinates")

plt.show()

# sol = lumi.load_sol(path)

# i = (len(frequencies) - 1) // 2
# φ, θ = np.meshgrid(
#     np.linspace(0, φmax, round(φmax / angres) + 1),
#     np.linspace(0, θmax, round(θmax / angres) + 1),
# )
# r = sol["dBi"]["o2"][
#     :, :, i
# ]  # or "o2@H" and "o2@V" for horizontal and vertical polarizations

# z = r * np.cos(θ)
# x = r * np.sin(θ) * np.cos(φ)
# y = r * np.sin(θ) * np.sin(φ)

# # Create the 3D plot
# fig = plt.figure(figsize=(10, 8))
# ax = fig.add_subplot(111, projection="3d")

# # Plot the surface
# ax.plot_surface(x, y, z, cmap="viridis")

# # Set labels
# ax.set_xlabel("X")
# ax.set_ylabel("Y")
# ax.set_zlabel("Z")
# ax.set_title("3D Surface Plot in Polar Coordinates")

# plt.show()

# fig = plt.figure()
# x = frequencies
# y = [10 * np.log10(v) for v in lumi.query(sol, "T1,1")]
# # y = [10 * math.log10(v) for v in lumi.query(sol, "To2@0,o1@0")]

# plt.plot(x, y)
# plt.xlabel("f [GHz]")
# plt.ylabel("S11 [dB]")
# plt.show()

# c.add_polygon(
#     [
#         (-l_patch / 2, -w_patch / 2),
#         (-l_patch / 2 + l_flap, -w_patch / 2),
#         (-l_patch / 2 + l_flap, -w_patch / 2 + w_flap),
#         (-l_patch / 2, -w_patch / 2 + w_flap),
#     ],
#     layer=PATCH,
# )
# c.add_polygon(
#     [
#         (-l_patch / 2, w_patch / 2),
#         (-l_patch / 2 + l_flap, w_patch / 2),
#         (-l_patch / 2 + l_flap, w_patch / 2 - w_flap),
#         (-l_patch / 2, w_patch / 2 - w_flap),
#     ],
#     layer=PATCH,
# )
# c.add_polygon(
#     [
#         (l_patch / 2, -w_patch / 2),
#         (-l_patch / 2 + l_flap, -w_patch / 2),
#         (-l_patch / 2 + l_flap, w_patch / 2),
#         (l_patch / 2, w_patch / 2),
#     ],
#     layer=PATCH,
# )

# c.add_polygon(
#     [
#         (-l_sub / 2, -w_sub / 2),
#         (-l_sub / 2, w_sub / 2),
#         (l_sub / 2, w_sub / 2),
#         (l_sub / 2, -w_sub / 2),
#     ],
#     layer=SUB,
# )

# line = c << gf.components.straight(
#     (l_sub - l_patch) / 2 + l_flap + margin, width=w_line
# )
# line.movex(-l_sub / 2 - margin)
# c.add_port(
#     f"o1", center=(-l_patch / 2 - l_line, 0), width=w_line, orientation=180, layer=WG
# )

# c.add_polygon(
#     [
#         (-(l_sub) / 2 - margin, -w_line / 2 - lateral_port_margin),
#         (-(l_sub) / 2 - margin, w_line / 2 + lateral_port_margin),
#         (-l_sub / 2, w_line / 2 + lateral_port_margin),
#         (-l_sub / 2, -w_line / 2 - lateral_port_margin),
#     ],
#     layer=SUB,
# )

# c << gf.components.bbox(
#     component=c, layer=BBOX, top=margin, bottom=margin, right=margin, left=0
# )
