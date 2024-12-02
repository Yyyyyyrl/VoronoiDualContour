#!/usr/bin/env python3

import sys
import re
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def parse_vd_info(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Initialize data structures
    voronoi_vertices = {}
    voronoi_facets = {}
    voronoi_cells = []
    isosurface_vertices = {}

    # Initialize current section
    current_section = None

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line == 'VoronoiVertices:':
            current_section = 'VoronoiVertices'
            i += 1
            continue
        elif line == 'VoronoiEdges:':
            current_section = 'VoronoiEdges'
            i += 1
            continue
        elif line == 'VoronoiVertexValues:':
            current_section = 'VoronoiVertexValues'
            i += 1
            continue
        elif line == 'VoronoiFacets:':
            current_section = 'VoronoiFacets'
            i += 1
            continue
        elif line == 'VoronoiCells:':
            current_section = 'VoronoiCells'
            i += 1
            continue
        elif line == 'IsosurfaceVertices:':
            current_section = 'IsosurfaceVertices'
            i += 1
            continue
        elif line.startswith('VoronoiDiagram:'):
            # Not needed
            i += 1
            continue

        # Now parse based on current_section
        if current_section == 'VoronoiVertices':
            # Parse VoronoiVertices
            if line.startswith('Index '):
                index = int(line.split('Index ')[1].strip(':'))
                # Next line should be 'VoronoiVertex:'
                i += 1
                line = lines[i].strip()
                if line != 'VoronoiVertex:':
                    # Error
                    pass
                # Next line should be 'Vertex: x y z'
                i += 1
                line = lines[i].strip()
                if line.startswith('Vertex:'):
                    coords = line.split('Vertex:')[1].strip()
                    x, y, z = map(float, coords.split())
                    voronoi_vertices[index] = (x, y, z)
                else:
                    pass
            i += 1
            continue
        elif current_section == 'VoronoiFacets':
            # Parse VoronoiFacets
            if line.startswith('Index '):
                index = int(line.split('Index ')[1].strip(':'))
                # Next line should be 'VoronoiFacet:'
                i += 1
                line = lines[i].strip()
                if line != 'VoronoiFacet:':
                    pass
                # Next lines should be:
                # Vertices indices: ...
                # Vertex values: ...
                i += 1
                line = lines[i].strip()
                if line.startswith('Vertices indices:'):
                    vertices_indices = list(map(int, line.split(':')[1].strip().split()))
                else:
                    pass
                i += 1
                line = lines[i].strip()
                if line.startswith('Vertex values:'):
                    vertex_values = list(map(float, line.split(':')[1].strip().split()))
                else:
                    pass
                voronoi_facets[index] = {'vertices_indices': vertices_indices, 'vertex_values': vertex_values}
            i += 1
            continue
        elif current_section == 'VoronoiCells':
            # Parse VoronoiCells
            if line == '':
                i += 1
                continue
            if line == 'VoronoiCell:':
                cell = {}
                # Read 'Cell index: i'
                i += 1
                line = lines[i].strip()
                if line.startswith('Cell index:'):
                    cell_index = int(line.split('Cell index:')[1].strip())
                    cell['cell_index'] = cell_index
                else:
                    pass
                # Read 'Delaunay vertex: x y z'
                i += 1
                line = lines[i].strip()
                if line.startswith('Delaunay vertex:'):
                    coords = line.split('Delaunay vertex:')[1].strip()
                    x, y, z = map(float, coords.split())
                    cell['delaunay_vertex'] = (x, y, z)
                else:
                    pass
                # Read 'Voronoi Vertices indices: ...'
                i +=1
                line = lines[i].strip()
                if line.startswith('Voronoi Vertices indices:'):
                    vertices_indices = list(map(int, line.split(':')[1].strip().split()))
                    cell['voronoi_vertices_indices'] = vertices_indices
                else:
                    pass
                # Read 'Facet indices: ...'
                i +=1
                line = lines[i].strip()
                if line.startswith('Facet indices:'):
                    facet_indices = list(map(int, line.split(':')[1].strip().split()))
                    cell['facet_indices'] = facet_indices
                else:
                    pass
                # Read 'IsoVertex start index: i'
                i +=1
                line = lines[i].strip()
                if line.startswith('IsoVertex start index:'):
                    isovertex_start_index = int(line.split(':')[1].strip())
                    cell['isovertex_start_index'] = isovertex_start_index
                else:
                    pass
                # Read 'Number of isoVertices: n'
                i +=1
                line = lines[i].strip()
                if line.startswith('Number of isoVertices:'):
                    number_of_isovertex = int(line.split(':')[1].strip())
                    cell['number_of_isovertex'] = number_of_isovertex
                else:
                    pass
                # Skip 'Cycles' section
                while i < len(lines):
                    line = lines[i].strip()
                    if line == '' or line.startswith('VoronoiCell:') or line.startswith('IsosurfaceVertices:'):
                        break
                    i += 1
                voronoi_cells.append(cell)
                continue
            i += 1
            continue
        elif current_section == 'IsosurfaceVertices':
            # Parse IsosurfaceVertices
            if line.startswith('Index '):
                index = int(line.split('Index ')[1].split(':')[0])
                coords = line.split(':')[1].strip()
                x, y, z = map(float, coords.split())
                isosurface_vertices[index] = (x, y, z)
            i += 1
            continue
        else:
            i +=1
            continue

    return voronoi_vertices, voronoi_facets, voronoi_cells, isosurface_vertices

def main():
    if len(sys.argv) < 3:
        print("Usage: python3 visVD.py vd_info.txt [cell indices]")
        print("Example: python3 visVD.py vd_info.txt [1 2 3]")
        sys.exit(1)

    filename = sys.argv[1]
    cell_indices = sys.argv[2:]
    # Remove brackets if present
    if cell_indices[0] == '[':
        cell_indices = cell_indices[1:]
    if cell_indices[-1] == ']':
        cell_indices = cell_indices[:-1]
    # Convert to integers
    cell_indices = [int(idx) for idx in cell_indices]

    # Parse the vd_info.txt file
    voronoi_vertices, voronoi_facets, voronoi_cells, isosurface_vertices = parse_vd_info(filename)

    # Create a 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for cell_index in cell_indices:
        # Find the cell with the given index
        cell = next((cell for cell in voronoi_cells if cell['cell_index'] == cell_index), None)
        if cell is None:
            print(f"Cell index {cell_index} not found.")
            continue

        # Get the Voronoi Vertices indices and coordinates
        vertex_indices = cell['voronoi_vertices_indices']
        verts = [voronoi_vertices[idx] for idx in vertex_indices]

        # Build a mapping from global vertex index to local index
        vertex_index_map = {global_idx: local_idx for local_idx, global_idx in enumerate(vertex_indices)}

        # Get the facets
        faces = []
        for facet_idx in cell['facet_indices']:
            if facet_idx not in voronoi_facets:
                continue
            facet = voronoi_facets[facet_idx]
            face_vertex_indices = facet['vertices_indices']
            # Map global indices to local indices
            face = []
            for vi in face_vertex_indices:
                if vi in vertex_index_map:
                    face.append(vertex_index_map[vi])
            if len(face) >= 3:
                face_coords = [verts[idx] for idx in face]
                faces.append(face_coords)

        # Plot the cell
        poly = Poly3DCollection(faces, alpha=0.5)
        # Optionally, set face color based on cell index
        poly.set_facecolor(np.random.rand(3,))
        ax.add_collection3d(poly)

        # Plot the vertices
        verts_array = np.array(verts)
        ax.scatter(verts_array[:,0], verts_array[:,1], verts_array[:,2], color='black', s=20)

        # Plot the iso vertex
        isovertex_index = cell['isovertex_start_index']
        if isovertex_index in isosurface_vertices:
            isovertex_coords = isosurface_vertices[isovertex_index]
            ax.scatter(isovertex_coords[0], isovertex_coords[1], isovertex_coords[2], color='red', s=100, marker='o')
        else:
            print(f"Isovertex index {isovertex_index} not found for cell {cell_index}.")

    # Set plot labels and show the plot
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.title('Voronoi Cells Visualization')
    plt.show()

if __name__ == "__main__":
    main()

