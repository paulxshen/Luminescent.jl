using Meshes
# 2D rectilinear grid
x = 0.0:0.2:1.0
y = [0.0, 0.1, 0.3, 0.7, 0.9, 1.0]
grid = RectilinearGrid(x, x, y)

n = Ngon((0, 0), (0.1, 0), (0.1, 0.1), (0, 0.1))
I = grid ∩ n