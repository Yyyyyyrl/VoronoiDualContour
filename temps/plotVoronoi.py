import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import csv

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

def plot_voronoi_3d(voronoi_vertices, voronoi_edges):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
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
    plt.show()

# Example usage
filename = 'VoronoiDiagram.csv'
voronoi_vertices, voronoi_edges = load_voronoi_from_csv(filename)
plot_voronoi_3d(voronoi_vertices, voronoi_edges)
