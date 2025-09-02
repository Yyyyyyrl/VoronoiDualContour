# Voronoi-Based Dual Contouring (VDC) - 3D Isosurface Extraction

## Overview

VDC is a C++ implementation of a Voronoi diagram-based dual contouring method for high-quality isosurface extraction from 3D scalar fields. The project leverages CGAL for computational geometry operations


## Installation

### Prerequisites
- CMake (version 3.10 or higher)
- CGAL (version 5.0 or higher)
- Boost libraries (for filesystem and program options)
- Python 3 (for visualization scripts only)

### Building from Source
```bash
git clone https://github.com/Yyyyyyrl/voronoi-dual-contouring.git
cd voronoi-dual-contouring
mkdir build && cd build
cmake ..
make
```

## Usage

### Basic Command Line Syntax
```bash
./vdc [OPTIONS] <isovalue> <(nhdr/nrrd) raw data file path>
```

### Full Options
```bash
./vdc -i <input> -o <output> -v <value> [options]

Required parameters:
  <file>          Input file path (NRRD format)
  <isovalue>      Isovalue for surface extraction

Processing options:
  -multi_isov     Enable multi-isosurface mode
  -single_isov    Enable separation of non-adjacent active cubes
  --supersample N Supersample input data by factor N
  --conv_H        Use convex hull method for Voronoi cells computation
  --sep_isov      Subsampling the input data for isosurface extraction

Output options:
  --out_csv <file>    Export Voronoi diagram to CSV file
  -off                Write the output in OFF format (default)
  -ply                Write the output in PLY format

```


## File Structure

## Detailed File Structure

```
.
├── CMakeLists.txt          # CMake build configuration
├── README.md               # Project documentation
├── compExec.py             # Output comparison script
├── visVD.py                # Mesh visualization tool
├── visV3d.py               # 3D visualization tool
├── include/                # Header files
│   ├── vdc.h               # Main header
│   ├── vdc_commandline.h   # Command line parsing
│   ├── vdc_delaunay.h      # Delaunay triangulation
│   ├── vdc_debug.h         # Debug utilities
│   ├── vdc_func.h          # Core algorithm functions
│   ├── vdc_grid.h          # Grid data structures
│   ├── vdc_io.h            # I/O operations
│   ├── vdc_type.h          # Type definitions
│   └── vdc_voronoi.h       # Voronoi diagram operations
└── src/                    # Source files
    ├── vdc.cpp             # Main application
    ├── vdc_commandline.cpp # Command line implementation
    ├── vdc_delaunay.cpp    # Delaunay implementation
    ├── vdc_debug.cpp       # Debug implementation
    ├── vdc_func_*.cpp      # Algorithm implementations
    ├── vdc_grid.cpp        # Grid operations
    ├── vdc_io.cpp          # I/O implementation
    └── vdc_voronoi.cpp     # Voronoi implementation
```

## API Documentation

The project provides comprehensive API documentation through Doxygen. To generate:
```bash
doxygen Doxyfile
```

