import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def read_grid_and_points_from_csv(filename):
    """
    Reads a CSV of the form:
      1) First line: nx, ny, nz, dx, dy, dz
      2) Second line: "x,y,z"
      3) Following lines: x,y,z
    Returns: (nx, ny, nz, dx, dy, dz, points_array)
    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    # Strip whitespace
    lines = [line.strip() for line in lines]

    # 1) Parse the first line for grid info
    grid_info_str = lines[0]  # e.g. "40,30,20,1.0,1.0,1.0"
    nx, ny, nz, dx, dy, dz = grid_info_str.split(',')
    nx, ny, nz = int(nx), int(ny), int(nz)
    dx, dy, dz = float(dx), float(dy), float(dz)

    # 2) The second line is likely "x,y,z" (header). We can skip it or parse it
    # We'll skip it:
    header_line = lines[1]

    # 3) The remaining lines are the point data
    data_lines = lines[2:]  # from the third line onward
    points_list = []
    for line in data_lines:
        x_str, y_str, z_str = line.split(',')
        x_val, y_val, z_val = float(x_str), float(y_str), float(z_str)
        points_list.append((x_val, y_val, z_val))

    # Convert to numpy array of shape (N, 3)
    points_array = np.array(points_list)

    return nx, ny, nz, dx, dy, dz, points_array

def generate_facets_as_lines(nx, ny, nz, dx, dy, dz):
    """
    Creates lines for the 6 bounding planes in 3D:
    x=0, x=(nx-1)*dx, y=0, y=(ny-1)*dy, z=0, z=(nz-1)*dz.
    
    """
    x_max = (nx - 1) * dx
    y_max = (ny - 1) * dy
    z_max = (nz - 1) * dz
    
    lines = []

    # We create the corners for each plane as a rectangle of 4 corners:
    # e.g. for x=0 plane: (0,0,0), (0,y_max,0), (0,y_max,z_max), (0,0,z_max)
    # Then close the loop.

    # x=0
    corners = np.array([
        [0, 0, 0],
        [0, y_max, 0],
        [0, y_max, z_max],
        [0, 0, z_max]
    ])
    lines.append(corners)

    # x=x_max
    corners = np.array([
        [x_max, 0, 0],
        [x_max, y_max, 0],
        [x_max, y_max, z_max],
        [x_max, 0, z_max]
    ])
    lines.append(corners)

    # y=0
    corners = np.array([
        [0, 0, 0],
        [x_max, 0, 0],
        [x_max, 0, z_max],
        [0, 0, z_max]
    ])
    lines.append(corners)

    # y=y_max
    corners = np.array([
        [0, y_max, 0],
        [x_max, y_max, 0],
        [x_max, y_max, z_max],
        [0, y_max, z_max]
    ])
    lines.append(corners)

    # z=0
    corners = np.array([
        [0, 0, 0],
        [x_max, 0, 0],
        [x_max, y_max, 0],
        [0, y_max, 0]
    ])
    lines.append(corners)

    # z=z_max
    corners = np.array([
        [0, 0, z_max],
        [x_max, 0, z_max],
        [x_max, y_max, z_max],
        [0, y_max, z_max]
    ])
    lines.append(corners)

    return lines

def main():

    if len(sys.argv) < 2:
        print("Usage: python(3) plotDP.py [filename].csv")
        sys.exit(1)

    # CSV file contains the dummy points and grid info
    filename = sys.argv[1]

    nx, ny, nz, dx, dy, dz, points_array = read_grid_and_points_from_csv(filename)
    xs = points_array[:, 0]
    ys = points_array[:, 1]
    zs = points_array[:, 2]

    # Create a 3D figure
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title("Dummy Points + Facets from" + filename)

    # Plot the dummy points in red
    ax.scatter(xs, ys, zs, c='red', s=10, label='Dummy Points')

    # Plot the bounding 6 facets as wireframe
    lines = generate_facets_as_lines(nx, ny, nz, dx, dy, dz)
    for corners in lines:
        # corners has shape (4,3), we add the first corner at the end to close loop
        closed_corners = np.vstack([corners, corners[0]])
        ax.plot(closed_corners[:, 0], closed_corners[:, 1], closed_corners[:, 2],
                color='blue', linewidth=1)

    # Label axes
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.legend()
    plt.show()

if __name__ == "__main__":
    main()
