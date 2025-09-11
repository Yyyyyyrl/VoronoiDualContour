import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import Voronoi

# Given points for the 3D Voronoi diagram
""" points = np.array([
    [0.5, 0.5, 0.5],
    [0.5, 0.5, 1.5],
    [0.5, 1.5, 0.5],
    [0.5, 1.5, 1.5],
    [1.5, 0.5, 0.5],
    [1.5, 0.5, 1.5],
    [1.5, 1.5, 0.5],
    [1.5, 1.5, 1.5]
]) """

points = np.array([
[0.5, 0.5, 0.5],
[0.5, 0.5, 1.5],
[0.5, 0.5, 2.5],
[0.5, 0.5, 3.5],
[0.5, 1.5, 0.5],
[0.5, 1.5, 1.5],
[0.5, 1.5, 2.5],
[0.5, 2.5, 0.5],
[0.5, 2.5, 1.5],
[0.5, 3.5, 0.5],
[1.5, 0.5, 0.5],
[1.5, 0.5, 1.5],
[1.5, 0.5, 2.5],
[1.5, 1.5, 0.5],
[1.5, 1.5, 1.5],
[1.5, 2.5, 0.5]
])

# Compute the Voronoi diagram
vor = Voronoi(points)

# Create a 3D plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot the points
ax.scatter(points[:, 0], points[:, 1], points[:, 2], c='b')

# Plot the Voronoi vertices
ax.scatter(vor.vertices[:, 0], vor.vertices[:, 1], vor.vertices[:, 2], c='r')

# Plot the Voronoi regions
for ridge in vor.ridge_vertices:
    if all(v >= 0 for v in ridge):
        ax.plot3D(vor.vertices[ridge, 0], vor.vertices[ridge, 1], vor.vertices[ridge, 2], 'g-')

# Counting the number of Voronoi vertices and edges
num_vertices = len(vor.vertices)
num_edges = sum(1 for simplex in vor.ridge_vertices if np.all(np.array(simplex) >= 0))

print("# vertices:" + str(num_vertices))
print("# edges: " + str(num_edges))


ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('3D Voronoi Diagram')

plt.show()
