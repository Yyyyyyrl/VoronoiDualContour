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
├── CMakeLists.txt            # CMake build configuration
├── README.md                 # Project documentation
├── cmake/                    # CMake helpers (if any)
├── include/                  # Header files
│   ├── core/                 # Core types, utilities, CLI
│   │   ├── vdc.h
│   │   ├── vdc_type.h
│   │   ├── vdc_utilities.h
│   │   ├── vdc_debug.h
│   │   └── vdc_commandline.h
│   ├── grid/                 # Grid data structures
│   │   └── vdc_grid.h
│   ├── io/                   # IO helpers (NRRD)
│   │   └── vdc_io.h
│   ├── delaunay/             # Delaunay types/utilities
│   │   └── vdc_delaunay.h
│   ├── voronoi/              # Voronoi structures and ops
│   │   └── vdc_voronoi.h
│   ├── algo/                 # High-level algorithm declarations
│   │   └── vdc_func.h
│   └── test/                 # Test helpers
│       └── test_vor.h
├── src/                      # Source files
│   ├── app/                  # Application entry points
│   │   └── vdc.cpp
│   ├── core/                 # Core utilities and CLI
│   │   ├── vdc_commandline.cpp
│   │   ├── vdc_debug.cpp
│   │   └── vdc_utilities.cpp
│   ├── grid/
│   │   └── vdc_grid.cpp
│   ├── io/
│   │   └── vdc_io.cpp
│   ├── algo/                 # Algorithm implementations
│   │   ├── vdc_func_delaunay.cpp
│   │   ├── vdc_func_helper.cpp
│   │   ├── vdc_func_isosurface.cpp
│   │   ├── vdc_func_modcyc.cpp
│   │   └── vdc_func_voronoi.cpp
│   ├── voronoi/              # Voronoi core and utilities
│   │   ├── vdc_voronoi_collapse.cpp
│   │   ├── vdc_voronoi_core.cpp
│   │   └── vdc_voronoi_matching.cpp
│   └── tests/                # Test drivers
│       ├── test_vor.cpp
│       └── test_modcyc.cpp
├── tools/
│   └── plot_modcyc.py        # Visualize modify-cycles test cases
├── data/                     # Sample volumes (NRRD/NHDR/RAW) and the outputs
├── compExec.py               # Output comparison helper
├── visVD.py                  # Voronoi diagram visualization
├── visV3d.py                 # 3D viewer helpers
├── visDP.py                  # Data-plane visualization
├── build/                    # Local build dir (generated)
└── buildLocal/               # Local build dir (generated)
```

## API Documentation

The project provides comprehensive API documentation through Doxygen. To generate:
```bash
doxygen Doxyfile
```
