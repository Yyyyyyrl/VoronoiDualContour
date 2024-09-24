import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import csv
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def load_voronoi_from_csv(filename):
    voronoi_vertices = []
    voronoi_edges = []
    with open(filename, newline='') as csvfile:
        reader = csv.reader(csvfile)
        section = None
        for row in reader:
            if row[0] == 'vertices':
                section = 'vertices'
            elif row[0] == 'edges':
                section = 'edges'
            elif section == 'vertices':
                # Parse vertices
                vertex = [float(coord) for coord in row]
                voronoi_vertices.append(vertex)
            elif section == 'edges':
                # Parse edges
                edge_type = row[0]
                if edge_type == 'Segment3' or edge_type == 'Line3':
                    p1 = [float(coord) for coord in row[1:4]]
                    p2 = [float(coord) for coord in row[4:7]]
                    voronoi_edges.append({'type': edge_type, 'points': [p1, p2]})
                elif edge_type == 'Ray3':
                    p1 = [float(coord) for coord in row[1:4]]
                    direction = [float(coord) for coord in row[4:7]]
                    voronoi_edges.append({'type': edge_type, 'points': [p1, direction]})
    return voronoi_vertices, voronoi_edges

def read_triangles_3d(filename):
    triangles = []
    vertices = {}
    with open(filename, 'r') as file:
        lines = file.readlines()
        current_poly = None
        for line in lines:
            line = line.strip()
            if line.startswith('Poly'):
                current_poly = []
                triangles.append(current_poly)
            elif line.startswith('Vertices:'):
                continue
            elif line:
                parts = line.split()
                vertex_id = int(parts[0])
                coords = tuple(float(x) for x in parts[1].strip('()').split(','))
                vertices[vertex_id] = coords
                current_poly.append(vertex_id)
    return vertices, triangles

def plot_triangles_3d(ax, vertices, triangles, color='red', label='Triangles'):
    for triangle in triangles:
        triangle_points = [vertices[vid] for vid in triangle]
        poly = [[triangle_points[0], triangle_points[1], triangle_points[2]]]
        ax.add_collection3d(Poly3DCollection(poly, facecolors=color, edgecolors='k', linewidths=1, alpha=0.5))

def plot_voronoi_3d(voronoi_vertices, voronoi_edges, ax):
    # Plot vertices
    for vertex in voronoi_vertices:
        ax.scatter(vertex[0], vertex[1], vertex[2], c='b', marker='o')

    # Plot edges
    for edge in voronoi_edges:
        if edge['type'] == 'Segment3':
            p1, p2 = edge['points']
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]], 'r-')
        elif edge['type'] == 'Line3':
            p1, p2 = edge['points']
            # Extend the line over a large interval
            direction = np.array(p2) - np.array(p1)
            line = np.array([p1 - 2 * direction, p2 + 2 * direction])
            ax.plot(line[:, 0], line[:, 1], line[:, 2], 'g-')
        elif edge['type'] == 'Ray3':
            p1, direction = edge['points']
            # Extend the ray from its origin
            line = np.array([p1, p1 + 2 * np.array(direction)])
            ax.plot(line[:, 0], line[:, 1], line[:, 2], 'y-')

    # Label axes
    ax.set_xlabel('X Axis')
    ax.set_ylabel('Y Axis')
    ax.set_zlabel('Z Axis')

def plot_combined_voronoi_and_triangles_3d(voronoi_file, triangles_file):
    voronoi_vertices, voronoi_edges = load_voronoi_from_csv(voronoi_file)
    vertices, triangles = read_triangles_3d(triangles_file)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Plot Voronoi diagram
    plot_voronoi_3d(voronoi_vertices, voronoi_edges, ax)

    # Plot additional triangles
    plot_triangles_3d(ax, vertices, triangles, color='red', label='Triangles')

    plt.show()

# Example usage
voronoi_file = input("CSV file path for Voronoi diagram: ")
triangles_file = input("txt file for vertices and triangles: ")  # Path to the triangles file
plot_combined_voronoi_and_triangles_3d(voronoi_file, triangles_file)
