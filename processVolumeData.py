import numpy as np
import nibabel as nib
import os
import sys

def parse_nhdr_header(filename):
    """
    Parse the NHDR header file to extract metadata about the scalar grid.

    Parameters:
        filename (str): The path to the NHDR file.

    Returns:
        dict: A dictionary containing the metadata from the NHDR file.
    """
    header = {}
    with open(filename, 'r') as file:
        for line in file:
            if ':' in line:
                key, value = line.strip().split(':')
                header[key.strip()] = value.strip()
    return header

def load_raw_data(header, raw_file_path):
    """
    Load the raw binary data from a file based on the NHDR metadata.

    Parameters:
        header (dict): Metadata containing data type and dimensions.
        raw_file_path (str): Path to the raw data file.

    Returns:
        np.array: A NumPy array containing the structured raw data.
    """
    dtype_map = {'unsigned char': np.uint8, 'short': np.int16, 'float': np.float32}
    data_type = dtype_map[header['type']]
    dimensions = tuple(map(int, header['sizes'].split()))
    with open(raw_file_path, 'rb') as file:
        data = np.fromfile(file, dtype=data_type).reshape(dimensions)
    return data

def load_nhdr_data(nhdr_path):
    """
    Load the raw binary data for NHDR files.

    Parameters:
        nhdr_path (str): Path to the NHDR file.

    Returns:
        np.array: A NumPy array containing the structured raw data.
    """
    header = parse_nhdr_header(nhdr_path)
    raw_file_path = os.path.join(os.path.dirname(nhdr_path), header['data file'])
    return load_raw_data(header, raw_file_path)


def read_nrrd_header(file_path):
    """Reads the header of an NRRD file to determine data dimensions and type."""
    header = {}
    with open(file_path, 'r') as file:
        for line in file:
            if line.strip() == '':
                break
            if ':' in line:
                key, value = line.strip().split(':', 1)
                header[key.strip()] = value.strip()
    return header

def read_nrrd_data(file_path):
    """Reads ASCII encoded NRRD data from a given file."""
    header = read_nrrd_header(file_path)
    dimensions = tuple(map(int, header['sizes'].split()))
    dtype_map = {
        'unsigned char': np.uint8,
        'short': np.int16,
        'float': np.float32
    }
    data_type = dtype_map[header['type']]
    
    data = []
    with open(file_path, 'r') as file:
        # Skip the header lines
        while True:
            line = file.readline()
            if line.strip() == '':  # Empty line indicates end of header
                break
        # Read the data
        for line in file:
            data.extend(map(float, line.split()))

    print(np.array(data).reshape(dimensions).astype(data_type))
    
    return np.array(data).reshape(dimensions).astype(data_type)
    
def load_data(file_path):
    if file_path.endswith('.nhdr'):
        return load_nhdr_data(file_path)
    elif file_path.endswith('.nrrd'):
        return read_nrrd_data(file_path)
    else:
        raise ValueError("Unsupported File Format")
    
def find_active_cubes(grid, isovalue):
    """
    Identify active cubes within the scalar grid that contain at least one bipolar edge.

    Parameters:
        grid (np.array): The 3D array of scalar values.
        isovalue (float): The threshold value to determine bipolar edges.

    Returns:
        list: A list of tuples representing the indices of active cubes.
    """
    nx, ny, nz = grid.shape
    active_cubes = []
    for i in range(nx - 1):
        for j in range(ny - 1):
            for k in range(nz - 1):
                if is_cube_active(grid, i, j, k, isovalue):
                    active_cubes.append((i, j, k))
    return active_cubes

def is_cube_active(grid, x, y, z, isovalue):
    """
    Determines if a cube has at least one bipolar edge.
    A bipolar edge is one where the scalar values at its endpoints straddle the isovalue.

    Parameters:
        grid (np.array): The 3D array of scalar values.
        x, y, z (int): Indices of the lower corner of the cube in the grid.
        isovalue (float): The threshold value for determining bipolar edges.

    Returns:
        bool: True if at least one edge is bipolar, False otherwise.
    """
    # Edges are defined by pairs of indices into the cube's corners
    # There are 12 edges in a cube
    edges = [
        ((x, y, z), (x + 1, y, z)),
        ((x, y, z), (x, y + 1, z)),
        ((x, y, z), (x, y, z + 1)),
        ((x + 1, y + 1, z), (x, y + 1, z)),
        ((x + 1, y + 1, z), (x + 1, y, z)),
        ((x + 1, y + 1, z), (x + 1, y + 1, z + 1)),
        ((x, y + 1, z + 1), (x, y, z + 1)),
        ((x, y + 1, z + 1), (x, y + 1, z)),
        ((x, y + 1, z + 1), (x + 1, y + 1, z + 1)),
        ((x + 1, y, z + 1), (x + 1, y, z)),
        ((x + 1, y, z + 1), (x, y, z + 1)),
        ((x + 1, y, z + 1), (x + 1, y + 1, z + 1))
    ]

    # Check each edge to determine if it is bipolar
    for (x0, y0, z0), (x1, y1, z1) in edges:
        val0 = grid[x0, y0, z0]
        val1 = grid[x1, y1, z1]
        if (val0 < isovalue and val1 > isovalue) or (val0 > isovalue and val1 < isovalue):
            return True  # This edge is bipolar

    return False  # No bipolar edges found


def process_volume_data(file_path, isovalue, output_file_path):
    """
    Process the volume data from NHDR or NRRD file, identify active cubes, compute centers, and write them to a file.

    Parameters:
        file_path (str): Path to the NHDR/NRRD file.
        isovalue (float): The isovalue to identify active cubes.
        output_file_path (str): Path to the output file where centers will be written.
    """
    grid = load_data(file_path)
    print("Raw data size:")
    print(grid.shape)
    
    active_cubes = find_active_cubes(grid, isovalue)
    centers = [(x + 0.5, y + 0.5, z + 0.5) for x, y, z in active_cubes]
    
    # Writing centers to file
    with open(output_file_path, "w") as file:
        for center in centers:
            file.write(f"{center[0]} {center[1]} {center[2]}\n")
    
    print("Number of active cubes:")
    print(len(centers))

def main():
    if len(sys.argv) < 3:
        print("Usage: python processVolumeData.py <file_path> <isovalue> <output_file_path>\n")
        sys.exit(1)
    
    file_path = sys.argv[1]
    isovalue = float(sys.argv[2])
    output_file_path = sys.argv[3]
    process_volume_data(file_path, isovalue, output_file_path)

if __name__ == "__main__":
    main()
